import math


# Taken from https://stackoverflow.com/questions/7237004/extending-a-line-segment-to-fit-into-a-bounding-box
def extend_line_to_bounds(point_1, point_2, bounds=(-1, -1, 1, 1)):

    # If we imagine extending the line until it crosses the top wall of the
    # bbox at point `(x_min, y_for_x_min)` and then imagine drawing
    # perpendicular lines from each point `(x_1, y_1)`, `(x_2, y_2)` to the wall
    # of the bbox, we end up with 2 perpendicular triangles with the same
    # angles - similar triangles. The rule of the similar triangles is that
    # the side lengths of two similar triangles are proportional.
    # That's how we get the equal ratios:
    # `| y_for_x_min - y_1 | / | x_min - x_1 | == | y_2 - y_1 | / | x_2 - x_1 |`
    # After we move some numbers from one to the other side of this equation,
    # we get the value for `y_for_x_min`. That's where the line should cross
    # the top wall of the bbox. We do the same for all other coordinates.
    # NOTE: These calculations are valid if one starts to draw a line from top
    # to bottom and from left to right. In case the direction is reverted, we
    # need to switch the min and max for each point (x, y). We do that below.

    # Extract
    x_1, y_1 = point_1[0], point_1[1]
    x_2, y_2 = point_2[0], point_2[1]
    x_min, y_min, x_max, y_max = bounds
    print(x_min)
    print(y_min)
    print(x_max)
    print(y_max)
    # Calculate
    y_for_x_min = y_1 + (y_2 - y_1) * (x_min - x_1) / (x_2 - x_1)
    y_for_x_max = y_1 + (y_2 - y_1) * (x_max - x_1) / (x_2 - x_1)
    x_for_y_min = x_1 + (x_2 - x_1) * (y_min - y_1) / (y_2 - y_1)
    x_for_y_max = x_1 + (x_2 - x_1) * (y_max - y_1) / (y_2 - y_1)

    # The line is vertical
    if (x_2 - x_1) < (y_2 - y_1):
        # The line is drawn from right to left
        if x_1 > x_2:
            # Switch the min and max x coordinates for y,
            # because the direction is from right (min) to left (max)
            y_for_x_min, y_for_x_max = y_for_x_max, y_for_x_min
    # The line is horizontal
    else:
        # The line is drawn from bottom to top
        if y_1 > y_2:
            # Switch the min and max y coordinates for x,
            # because the direction is from bottom (min) to top (max)
            x_for_y_min, x_for_y_max = x_for_y_max, x_for_y_min

    # The line is drawn from right to left
    if x_1 > x_2:
        # Get the maximal value for x_1.
        # When `x_for_y_min < x_min`(line goes out of the
        # bbox from the top), we clamp to x_min.
        x_1 = max(max(int(x_for_y_min), x_min), x_1)
    # The line is drawn from left to right
    else:
        # Get the minimal value for x_1.
        # When `x_for_y_min < x_min`(line goes out of the
        # bbox from the top), we clamp to x_min.
        if math.isinf(x_for_y_min):
            x_1 = min(x_min, x_1)
        else:
            x_1 = min(max(int(x_for_y_min), x_min), x_1)

    # Get the maximal value for x_2.
    # When `x_for_y_max > x_max` (line goes out of the
    # bbox from the bottom), we clamp to x_max.
    if math.isinf(x_for_y_max):
        x_2 = max(x_max, x_2)
    else:
        x_2 = max(min(int(x_for_y_max), x_max), x_2)

    # Get the minimal value for y_1
    # When `y_for_x_min < y_min`(line goes out of the
    # bbox from the left), we clamp to y_min.
    if math.isinf(y_for_x_min):
        y_1 = min(y_min, y_max)
    else:
        y_1 = min(max(int(y_for_x_min), y_min), y_max)

    # Get the minimal value for y_2
    if math.isinf(y_for_x_min):
        y_2 = y_max
    else:
        y_2 = min(int(y_for_x_max), y_max)
    # Done
    return [x_1, y_1, x_2, y_2]


def extend_line(xmin, ymin, xmax, ymax, x1, y1, x2, y2):

    """
    Extend a line so that it reaches the walls of the bbox.

    Args:
        xmin(int): The very left coordinate of the bbox.
        ymin(int): The very top coordinate of the bbox.
        xmax(int): The very right coordinate of the bbox.
        ymax(int): The very bottom coordinate of the bbox.
        x1(int): The start x coordinate of the line.
        y1(int): The start y coordinate of the line.
        x2(int): The end x coordinate of the line.
        y2(int): The end y coordinate of the line.

    Returns:
        - (list): The start and end (x, y) coordinates of the extended line.
    """

    # If we imagine extending the line until it crosses the top wall of the
    # bbox at point `(xmin, y_for_xmin)` and then imagine drawing
    # perpendicular lines from each point `(x1, y1)`, `(x2, y2)` to the wall
    # of the bbox, we end up with 2 perpendicular trianlges with the same
    # angles - similar triangles. The rule of the similar triangles is that
    # the side lengths of two similar triangles are proportional.
    # That's how we get the equal ratios:
    # `| y_for_xmin - y1 | / | xmin - x1 | == | y2 - y1 | / | x2 - x1 |`
    # After we move some numbers from one to the other side of this equation,
    # we get the value for `y_for_xmin`. That's where the line should cross
    # the top wall of the bbox. We do the same for all other coordinates.
    # NOTE: These calculations are valid if one starts to draw a line from top
    # to botton and from left to right. In case the direction is reverted, we
    # need to switch the min and max for each point (x, y). We do that below.
    y_for_xmin = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1)
    y_for_xmax = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1)
    x_for_ymin = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1)
    x_for_ymax = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1)

    # The line is vertical
    if (x2 - x1) < (y2 - y1):
        # The line is drawn from right to left
        if x1 > x2:
            # Switch the min and max x coordinates for y,
            # because the direction is from right (min) to left (max)
            y_for_xmin, y_for_xmax = y_for_xmax, y_for_xmin
    # The line is horizontal
    else:
        # The line is drawn from bottom to top
        if y1 > y2:
            # Switch the min and max y coordinates for x,
            # because the direction is from bottom (min) to top (max)
            x_for_ymin, x_for_ymax = x_for_ymax, x_for_ymin

    # The line is drawn from right to left
    if x1 > x2:
        # Get the maximal value for x1.
        # When `x_for_ymin < xmin`(line goes out of the
        # bbox from the top), we clamp to xmin.
        x1 = max(max(int(x_for_ymin), xmin), x1)
    # The line is drawn from left to right
    else:
        # Get the minimal value for x1.
        # When `x_for_ymin < xmin`(line goes out of the
        # bbox from the top), we clamp to xmin.
        if math.isinf(x_for_ymin):
            x1 = min(xmin, x1)
        else:
            x1 = min(max(int(x_for_ymin), xmin), x1)

    # Get the maximal value for x2.
    # When `x_for_ymax > xmax` (line goes out of the
    # bbox from the bottom), we clamp to xmax.
    if math.isinf(x_for_ymax):
        x2 = max(xmax, x2)
    else:
        x2 = max(min(int(x_for_ymax), xmax), x2)

    # Get the minimal value for y1
    # When `y_for_xmin < ymin`(line goes out of the
    # bbox from the left), we clamp to ymin.
    if math.isinf(y_for_xmin):
        y1 = min(ymin, ymax)
    else:
        y1 = min(max(int(y_for_xmin), ymin), ymax)

    # Get the minimal value for y2
    if math.isinf(y_for_xmin):
        y2 = ymax
    else:
        y2 = min(int(y_for_xmax), ymax)
    # Done
    return [x1, y1, x2, y2]