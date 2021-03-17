import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import math as m
from descartes import PolygonPatch
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd

def newPoint(lat, lon, direction, distance):
    '''
    Inputs:  initial lat, lon pair, a distance (feet), and a direction of travel
    Outputs: new lat, lon pair after translating `distance` in `direction`

    Inspired from: https://www.movable-type.co.uk/scripts/latlong.html
    '''

    R = 6373.0  # Radius of the Earth
    brng = m.radians(direction)  # Bearing is 90 degrees, converted to radians
    d = distance / 3280.84  # Convert distance from feet to km
    lat1 = m.radians(lat)  # Current lat point converted to radians
    lon1 = m.radians(lon)  # Current lon point converted to radians
    lat2 = m.asin(m.sin(lat1) * m.cos(d / R) + m.cos(lat1) * m.sin(d / R) * m.cos(brng))
    lon2 = lon1 + m.atan2(m.sin(brng) * m.sin(d / R) * m.cos(lat1), m.cos(d / R) - m.sin(lat1) * m.sin(lat2))

    # convert to degrees
    lat2 = m.degrees(lat2)
    lon2 = m.degrees(lon2)

    return lat2, lon2

def distance(lat1, lon1, lat2, lon2):
    '''
    Inputs: 2 lat, lon pairs representing 2 points
    Outputs: distance (feet) between 2 points

    Applies Haversine formula
    '''

    R = 6373.0 # Radius of the Earth
    lat1 = m.radians(lat1)
    lon1 = m.radians(lon1)
    lat2 = m.radians(lat2)
    lon2 = m.radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = m.sin(dlat / 2) ** 2 + m.cos(lat1) * m.cos(lat2) * m.sin(dlon / 2) ** 2
    c = 2 * m.atan2(m.sqrt(a), m.sqrt(1 - a))
    distance = R * c * 3280.84  # returns distance in feet

    return distance

def angle(x1, y1, x2, y2):
    '''
    Inputs: 2 vectors
    Outputs: Angle between 2 vectors, [0, pi]

    Uses dot product to find angle between 2 vectors
    '''

    numer = (x1 * x2 + y1 * y2)
    denom = m.sqrt((x1 ** 2 + y1 ** 2) * (x2 ** 2 + y2 ** 2))
    numer = 1e-17 if numer == 0 else numer
    denom = 2e-17 if denom == 0 else denom
    if abs(numer / denom) < 1.00000001 and abs(numer / denom) > .99999999:
        denom += .00000001 # account for when numer == denom

    return m.acos(numer / denom) / m.pi * 180

def visualizeSinglePoly(polygon, proposal_points, plotProposals=True):
    '''
    Plots polygon and line proposals
    '''

    fig = plt.figure()
    ax = fig.gca()
    ax.add_patch(PolygonPatch(polygon, alpha=.5))
    plt.xticks(rotation=45)
    plt.title('Line Proposals')
    plt.xlabel('Longitude', labelpad=20)
    plt.ylabel('Latitude', labelpad=20)
    ax.axis('scaled')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.tight_layout()
    plt.grid(True)

    if plotProposals:
        for p1, p2 in proposal_points:
            lon_values = (p1[1], p2[1])
            lat_values = (p1[0], p2[0])
            plt.plot(lon_values, lat_values, '-', color='navy')

    plt.show()

    return

def getCoordsFromPoly(polygon):
    '''
    Extracts coordinates of points on Polygon
    '''

    coords = list(polygon.exterior.coords)
    coords = coords[:-1] # remove last point due to duplicate start/end point

    return coords

def getCorners(polygon):
    '''
    Extracts critical corners from polygon
    '''

    coords = getCoordsFromPoly(polygon)  # coords of polygon

    corners = []  # list of corners of polygon
    for ind, point in enumerate(coords):
        # for every point within perimeter of polygon, point[i], point[i-1], point[i-2]
        # are used to determine the angle and total distance between the 3 points

        p1 = coords[ind]
        ref = coords[ind - 1]
        p2 = coords[ind - 2]
        x1, y1 = p1[0] - ref[0], p1[1] - ref[1]
        x2, y2 = p2[0] - ref[0], p2[1] - ref[1]
        intermediate_angle = angle(x1, y1, x2, y2)
        intermediate_distance = m.hypot(p2[0] - ref[0], p2[1] - ref[1]) + m.hypot(ref[0] - p1[0], ref[1] - p1[1])

        if (intermediate_angle > 60 and intermediate_angle < 120):
            corners.append(ref)
        elif (intermediate_angle > 10 and intermediate_angle < 150 and intermediate_distance > 0.01 * polygon.length):
            corners.append(ref)

    return corners

def defineCriticalEdge(corners, orientation):
    '''
    Extracts critical edge of polygon
    Critical edge of polygon = side of polygon wrt which line proposals will be parallel
    Finds the slope and intercept in cartesian plane corresponding to the line that represents the
    critical edge of polygon that has an angle close to 90 degrees;
    This is used to determine the orientation (slope and intercept) of linestring placements
    p1 and p2 correspond to the 2 corners of the polygon that define the longest side
    '''

    longest_length = 0

    for ind, corner in enumerate(corners):
        if corners[ind][0] == corners[ind - 1][0]:  # two points have same x value, cannot be accurately fit using polyfit
            continue

        if m.hypot(corners[ind][1] - corners[ind - 1][1], corners[ind][0] - corners[ind - 1][0]) < .05 * polygon.length:
            # ensures proposed side of polygon is not < 5% of total perimeter of polygon
            continue

        m_iter, b_iter = np.polyfit([corners[ind - 1][0], corners[ind][0]], [corners[ind - 1][1], corners[ind][1]], 1)

        if m.hypot(corners[ind][0] - corners[ind - 1][0], corners[ind][1] - corners[ind - 1][1]) > longest_length:

            if abs((90 - np.arctan(m_iter) * 180 / m.pi) - orientation) < 50 or abs((90 - np.arctan(m_iter) * 180 / m.pi) - (orientation + 180)) < 50:

                longest_length = m.hypot(corners[ind][0] - corners[ind - 1][0], corners[ind][1] - corners[ind - 1][1])
                m_longest, b_longest = m_iter, b_iter
                p1 = corners[ind - 1]
                p2 = corners[ind]

    return m_longest, b_longest, p1, p2


def generateCriticalEdgePoints(m_longest, b_longest):
    '''
    Generate points along critical edge (x, y)
    inv_slop_angle: finds the angle perpendicular to critical edge; the linestock proposals will be offset in this direction
    '''

    if abs(90 - np.arctan(m_longest) * 180 / m.pi - 0) < .1 or abs(90 - np.arctan(m_longest) * 180 / m.pi - 180) < .1:
        resolution = .00000005
    elif abs(90 - np.arctan(m_longest) * 180 / m.pi - 0) < 10 or abs(90 - np.arctan(m_longest) * 180 / m.pi - 180) < 10:
        resolution = .000001
    else:
        resolution = .00001

    x = np.arange(polygon.bounds[0], polygon.bounds[2], resolution)
    y = m_longest * x + b_longest
    inv_slope_angle = 90 - np.arctan(-1 / m_longest) * 180 / m.pi

    return x, y, inv_slope_angle



def linestringProposal(polygon, orientation, spacing, min_edge_spacing = 330, min_len = 500, plots_on = False, adjustment_factor = 0):
    '''
    polygon: a single shapely.geometry.Polygon object
    orientation: general orientation that proposed linestring should lay, in degrees [0, 360)
    spacing: distance between proposed linestring, feet
    min_edge_spacing: distance offset between any proposed linestring at polygon edge
    min_len: minimum length of proposed linestring, feet
    plots_on: for visualizations
    adjustment_factor: for use when centering linestring proposals within polygon, or
        for irregular-shaped polygons with custom linestring spacing requirements

    returns: proposals, a list of endpoints of proposed linestrings
    '''

    corners = getCorners(polygon)

    if plots_on: # plots corners
        fig = plt.figure()
        ax = fig.gca()
        ax.add_patch(PolygonPatch(polygon, alpha=.5))
        ax.axis('scaled')
        plt.title('Major Corners')
        plt.xlabel('Longitude', labelpad=20)
        plt.ylabel('Latitude', labelpad=20)
        plt.xticks(rotation=45)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        [plt.plot(corner[0], corner[1], 'ro') for corner in corners]
        plt.tight_layout()
        plt.grid(True)
        plt.show()

    m_longest, b_longest, p1, p2 = defineCriticalEdge(corners, orientation)
    x, y, inv_slope_angle = generateCriticalEdgePoints(m_longest, b_longest)

    # 0 represents p1
    # 0.5 represents midpoint between p1 and p2
    # 1 represents p2
    percentage_along_line = [-2, -1.75, -1.5, -1.25, -1, -.75, -.5, -.25, 0, .25, .5, .75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]
    p = {percentage: ((1-percentage) * p1[0] + percentage * p2[0], (1-percentage) * p1[1] + percentage * p2[1]) for percentage in percentage_along_line}
    edge_offset = max(spacing / 2, min_edge_spacing)  # distance linestrings must maintain away from polygon boundary

    # determines which side relative to edge of polygon linestrings will be proposed and finds midpoint of the first linestring
    lat_pos, lon_pos = newPoint(p[.5][1], p[.5][0], inv_slope_angle, edge_offset)
    lat_neg, lon_neg = newPoint(p[.5][1], p[.5][0], inv_slope_angle + 180, edge_offset)
    polygon_centroid_lat, polygon_centroid_lon = polygon.centroid.coords.xy[1][0], polygon.centroid.coords.xy[0][0]

    if distance(lat_pos, lon_pos, polygon_centroid_lat, polygon_centroid_lon) > distance(lat_neg, lon_neg, polygon_centroid_lat, polygon_centroid_lon):
        inv_slope_angle += 180

    lat_forward = {percentage: newPoint(p[percentage][1], p[percentage][0], inv_slope_angle, edge_offset)[0] for percentage in percentage_along_line}
    lon_forward = {percentage: newPoint(p[percentage][1], p[percentage][0], inv_slope_angle, edge_offset)[1] for percentage in percentage_along_line}
    lat_backward = {percentage: newPoint(lat_forward[percentage], lon_forward[percentage], inv_slope_angle + 180, spacing)[0] for percentage in percentage_along_line}
    lon_backward = {percentage: newPoint(lat_forward[percentage], lon_forward[percentage], inv_slope_angle + 180, spacing)[1] for percentage in percentage_along_line}

    edge_spacing_mod_1 = newPoint(lat_forward[.5], lon_forward[.5], inv_slope_angle, edge_offset)  # modified lat, lon to account for irregularities in polygon geometry
    edge_spacing_latlon_1 = m.hypot(edge_spacing_mod_1[0] - lat_forward[.5], edge_spacing_mod_1[1] - lon_forward[.5])
    edge_spacing_mod_2 = newPoint(lat_forward[.5], lon_forward[.5], inv_slope_angle + 90, edge_offset)  # modified lat, lon to account for irregularities in polygon geometry
    edge_spacing_latlon_2 = m.hypot(edge_spacing_mod_2[0] - lat_forward[.5], edge_spacing_mod_2[1] - lon_forward[.5])
    edge_spacing_latlon = np.mean([edge_spacing_latlon_1, edge_spacing_latlon_2])
    buffered = polygon.buffer(-.8 * edge_spacing_latlon)  # internally buffered geometry of polygon with all sides buffered by lease_edge distance

    if plots_on:  # plot internal buffer
        fig = plt.figure()
        ax = fig.gca()
        ax.add_patch(PolygonPatch(polygon, alpha=.5))
        ax.add_patch(PolygonPatch(buffered, alpha=.5))
        ax.axis('scaled')
        plt.title('Internal Buffer')
        plt.xlabel('Longitude', labelpad=20)
        plt.ylabel('Latitude', labelpad=20)
        plt.xticks(rotation=45)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        [plt.plot(corner[0], corner[1], 'ro') for corner in corners]
        plt.tight_layout()
        plt.grid(True)
        plt.show()

    lat_lon_bool_forward = [Point(lon_forward[percentage], lat_forward[percentage]).within(buffered) for percentage in percentage_along_line]
    lat_lon_bool_backward = [Point(lon_backward[percentage], lat_backward[percentage]).within(buffered) for percentage in percentage_along_line]
    proposals = []  # list of lengths of proposed linestrings and list of points on linestrings

    def iterateLinestrings(lat_lon_bool_list, lat_dict, lon_dict, forward):
        offset = 0 if forward else 180
        while True in lat_lon_bool_list:

            b_well = lat_dict[.5] - m_longest * lon_dict[.5]
            y_well = m_longest * x + b_well

            coords_all_linestring = []  # [[coords1], [corods2]] or [[(x, y)]]
            coords_linestring = []
            count_along_linestring = 0
            for ind, point in enumerate(x):
                if Point(x[ind], y_well[ind]).within(buffered):
                    coords_linestring.append((x[ind], y_well[ind]))
                if len(coords_linestring) != 0:
                    count_along_linestring += 1
                if count_along_linestring != len(coords_linestring):
                    coords_all_linestring.append(coords_linestring)
                    coords_linestring = []
                    count_along_linestring = 0

            for linestring in coords_all_linestring:
                proposals.append(linestring)

            lat_dict = {percentage: newPoint(lat_dict[percentage], lon_dict[percentage], inv_slope_angle + offset, spacing)[0] for percentage in percentage_along_line}
            lon_dict = {percentage: newPoint(lat_dict[percentage], lon_dict[percentage], inv_slope_angle + offset, spacing)[1] for percentage in percentage_along_line}
            lat_lon_bool_list = [Point(lon_dict[percentage], lat_dict[percentage]).within(buffered) for percentage in percentage_along_line]

        return lat_dict, lon_dict

    lat_backward, lon_backward = iterateLinestrings(lat_lon_bool_backward, lat_backward, lon_backward, forward=False)
    lat_forward, lon_forward = iterateLinestrings(lat_lon_bool_forward, lat_forward, lon_forward, forward=True)
    lat_forward = {percentage: newPoint(lat_forward[percentage], lon_forward[percentage], inv_slope_angle + 180, spacing)[0] for percentage in percentage_along_line}
    lon_forward = {percentage: newPoint(lat_forward[percentage], lon_forward[percentage], inv_slope_angle + 180, spacing)[1] for percentage in percentage_along_line}
    lateral_len_list = [distance(linestring[0][1], linestring[0][0], linestring[-1][1], linestring[-1][0]) for linestring in proposals]
    remove_lateral_ind = [ind for ind, linestring in enumerate(lateral_len_list) if linestring < min_len]

    # min_len check
    proposals = [linestring for ind, linestring in enumerate(proposals) if ind not in remove_lateral_ind]
    lateral_len_list = [linestring for ind, linestring in enumerate(lateral_len_list) if ind not in remove_lateral_ind]

    # Consolidate linestring to start and ends coords
    proposals = [((linestring[0][1], linestring[0][0]), (linestring[-1][1], linestring[-1][0])) for linestring in proposals]

    # Center linestring proposals within polygon
    center_distance = spacing / 50
    center_count = 0

    while Point(lon_forward[.5], lat_forward[.5]).within(polygon):
        lat_forward = {percentage: newPoint(lat_forward[percentage], lon_forward[percentage], inv_slope_angle, center_distance)[0] for
               percentage in percentage_along_line}
        lon_forward = {percentage: newPoint(lat_forward[percentage], lon_forward[percentage], inv_slope_angle, center_distance)[1] for
               percentage in percentage_along_line}
        center_count += 1

    distance_to_edge = (center_count - 1) * center_distance
    distance_to_shift = (distance_to_edge + edge_offset) / 2 - edge_offset

    proposals = [newPoint(coord[0], coord[1], inv_slope_angle, adjustment_factor * distance_to_shift) for linestring in proposals for coord in linestring]
    proposals = [(proposals[ind], proposals[ind + 1]) for ind in range(0, len(proposals), 2)]

    return proposals

if __name__ == "__main__":

    path = r'/Users/mwang/Documents/Active Projects/demo.shp'
    gdf = gpd.read_file(path)
    gdf = gdf.set_crs(epsg=4326) # set CRS

    for polygon_id in gdf.UNIT_NAME:

        print(polygon_id)
        polygon = gdf.loc[gdf.UNIT_NAME == polygon_id, 'geometry'].item()
        visualizeSinglePoly(polygon, [], plotProposals=False) # plot base polygon
        proposals = linestringProposal(polygon, orientation=90, spacing=500, min_len=800, min_edge_spacing=330, plots_on=True)
        visualizeSinglePoly(polygon, proposals, plotProposals=True) # plot polygon + proposals

