import csv
import sys
##### this performs the ray casting algorithm -  reads the convex hull mannered coordinates in a csv file and test points whidh
# are to be tested whether they are in the convex hull and stores as a output file  with third column as boolean whether it is inside or outside the convex hull
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def read_points_from_csv(filename):
    points = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 2:
                x = float(row[0])
                y = float(row[1])
                points.append(Point(x, y))

    print("length",len(points))
    return points

def is_inside_polygon(polygon, point):
    crossings = 0
    n = len(polygon)

    for i in range(n):
        p1 = polygon[i]
        p2 = polygon[(i + 1) % n]

        if ((p1.y <= point.y < p2.y) or (p2.y <= point.y < p1.y)) and (point.x < (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y) + p1.x):
            crossings += 1

    return crossings % 2 != 0

def write_output_to_csv(filename, test_points, convex_hull):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        for point in test_points:
            inside = is_inside_polygon(convex_hull, point)
            writer.writerow([point.x, point.y, int(inside)])

def main():
    convex_hull_file = sys.argv[1]
    test_points_file = 'nodes_coordinates.csv'
    output_file = 'nodes_marker.csv'

    convex_hull = read_points_from_csv(convex_hull_file)
    if len(convex_hull) < 3:
        print("At least 3 points are required to define the convex hull.")
        return

    test_points = read_points_from_csv(test_points_file)
    if len(test_points) == 0:
        print("No points found in the test points CSV file.")
        return

    write_output_to_csv(output_file, test_points, convex_hull)

if __name__ == '__main__':
    main()
