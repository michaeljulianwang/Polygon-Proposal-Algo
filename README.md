# Polygon-Proposal-Algo

Python3 algorithm that proposes linestring locations within an arbitrary-shaped polygon. 

Proposals account for user-specified constraints, including:
- `orientation`: desired cardinal
- `min_len`: minimum linestring length, below which proposals are dropped
- `spacing`: minimum distance between adjacent proposed linestrings
- `min_edge_spacing`: minimum distance proposed linestrings must maintain away from polygon exterior

Beyond above constraints, algorithm maximizes linestring count and length.

Starter demo data (shapefile) included for convenience.

Note: `adjustment_factor` needs to be used when centering linestrings within polygon, or when handling irregular-shaped polygons with custom linestring spacing requirements.

Note: demo data and images are arbitrarily wrt latitude/longitude.

# Algorithm
0. Pre-processing of polygons
1. Extract critical corners from polygon
2. Apply reverse buffer, respecting `min_edge_spacing`
3. Propose linestrings within buffered polygon wrt critical edge, respecting `spacing`

# Example 1 (polygon ID 21)

|  <img src="./img/21_polygon.png" alt="Polygon 21 Raw" width="500"> | 
|:--:| 
| Raw Polygon |

|  <img src="./img/21_corners.png" alt="Polygon 21 Corners" width="500"> | 
|:--:| 
| Critical Corners |

|  <img src="./img/21_buffer.png" alt="Polygon 21 Buffer" width="500"> | 
|:--:| 
| Internal Buffer |

|  <img src="./img/21_proposals_0.png" alt="Polygon 21 Proposals 0" width="500"> | 
|:--:| 
| Proposals, Orientation = 0&deg; |

|  <img src="./img/21_proposals_90.png" alt="Polygon 21 Proposals 90" width="500"> | 
|:--:| 
| Proposals, Orientation = 90&deg; |




# Example 2 (polygon ID 31)

|  <img src="./img/31_polygon.png" alt="Polygon 31 Raw" width="500"> | 
|:--:| 
| Raw Polygon |

|  <img src="./img/31_corners.png" alt="Polygon 31 Corners" width="500"> | 
|:--:| 
| Critical Corners |

|  <img src="./img/31_buffer.png" alt="Polygon 31 Buffer" width="500"> | 
|:--:| 
| Internal Buffer |

|  <img src="./img/31_proposals_135.png" alt="Polygon 31 Proposals 135" width="500"> | 
|:--:| 
| Proposals, Orientation = 135&deg; |

|  <img src="./img/31_proposals_45_a.png" alt="Polygon 31 Proposals 45_a" width="500"> | 
|:--:| 
| Proposals, Orientation = 45&deg;, Spacing = 300, MinLen = 100 |

|  <img src="./img/31_proposals_45_b.png" alt="Polygon 31 Proposals 45_b" width="500"> | 
|:--:| 
| Proposals, Orientation = 45&deg;, Spacing = 200, MinLen = 0 <br /> <em> (note tighter spacing density and inclusion of <br /> shorter linestring proposals, as compared to above) </em> |

