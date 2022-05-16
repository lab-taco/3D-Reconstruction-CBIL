from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.geometry import Polygon
from matplotlib import pyplot as plt

def point_dist_to_convex_hull(dist, plotting=False, plot_label='DIST.'):
    convex_hull = MultiPoint(dist).convex_hull
    plt.plot(*convex_hull.exterior.xy, label=plot_label)
    return convex_hull