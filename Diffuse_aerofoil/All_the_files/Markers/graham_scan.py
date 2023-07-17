import csv
# this program is used to convert the coordinates into anticlockwise manner of the convex hull and it should be fed into the markerprogram.py
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def orientation(p, q, r):
    val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    if val == 0:
        return 0  # Collinear
    return 1 if val > 0 else 2  # Clockwise or Counterclockwise

def compare_points(p0, p1, p2):
    o = orientation(p0, p1, p2)
    if o == 0:
        return (p0.x * p0.x + p0.y * p0.y) < (p1.x * p1.x + p1.y * p1.y)
    return o == 2  # Sort in counterclockwise order

def find_bottommost_point(points):
    bottommost = points[0]
    for p in points:
        if p.y < bottommost.y or (p.y == bottommost.y and p.x < bottommost.x):
            bottommost = p
    return bottommost

def transform_to_counterclockwise(points):
    bottommost = find_bottommost_point(points)
    modified_points = points.copy()
    modified_points.remove(bottommost)
    modified_points.sort(key=lambda p: (compare_points(bottommost, p, p), -(p.y - bottommost.y) / (p.x - bottommost.x) if p.x != bottommost.x else float('inf')))
    modified_points.insert(0, bottommost)
    return modified_points

def read_points_from_csv(filename):
    points = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) == 2:
                x, y = map(float, row)
                points.append(Point(x, y))
    return points

def write_points_to_csv(points, filename):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        for point in points:
            writer.writerow([point.x, point.y])
    print(f"Points written to file: {filename}")

def main():
    input_file = './points_not_in_convex_hull_manner/0012_25deg.csv'
    output_file = '0012_25deg_convex_hull.csv'
    points = read_points_from_csv(input_file)
    if not points:
        print("No points found in the CSV file!")
        return
    counterclockwise_points = transform_to_counterclockwise(points)
    write_points_to_csv(counterclockwise_points, output_file)

if __name__ == '__main__':
    main()
