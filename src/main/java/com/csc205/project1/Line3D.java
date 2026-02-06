package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a line segment in three-dimensional Euclidean space.
 *
 * This class demonstrates several fundamental object-oriented design patterns
 * and principles:
 *
 * 1. COMPOSITION PATTERN:
 *    - Composed of two Point3D objects (start and end points)
 *    - Delegates coordinate operations to Point3D class
 *    - "Has-a" relationship rather than inheritance
 *    - Demonstrates code reuse through composition over inheritance
 *
 * 2. VALUE OBJECT PATTERN:
 *    - Immutable design: All fields are final and private
 *    - No setters provided - transformations return new instances
 *    - Thread-safe by design due to immutability
 *    - Equality based on endpoints rather than identity
 *
 * 3. ENCAPSULATION:
 *    - Private fields with controlled access through public methods
 *    - Internal representation (two points vs. point + vector) is hidden
 *    - Clients interact through well-defined interface
 *
 * 4. SINGLE RESPONSIBILITY PRINCIPLE:
 *    - Focused solely on representing and manipulating 3D line segments
 *    - Geometric calculations specific to lines
 *    - Does not handle rendering, persistence, or other concerns
 *
 * 5. ALGORITHM FOUNDATIONS:
 *    - Parametric line equations
 *    - Point-line distance calculations (projection algorithms)
 *    - Line-line distance (closest point of approach)
 *    - Intersection detection (computational geometry)
 *    - Vector projections and rejections
 *
 * 6. DATA STRUCTURE PRINCIPLES:
 *    - Efficient representation using two points
 *    - O(1) access to endpoints
 *    - Immutable structure prevents concurrent modification issues
 *    - Lazy calculation of derived properties (direction vector)
 *
 * Mathematical Representation:
 * A line can be represented parametrically as: P(t) = start + t * direction
 * where t ∈ [0, 1] for a line segment, or t ∈ ℝ for an infinite line.
 *
 * @author Generated for Spring Framework style demonstration
 * @version 1.0
 */
public class Line3D {

    private static final Logger logger = Logger.getLogger(Line3D.class.getName());

    // Endpoints of the line segment
    private final Point3D start;
    private final Point3D end;

    // Tolerance for floating-point comparisons
    private static final double EPSILON = 1e-10;

    /**
     * Constructs a new Line3D with the specified start and end points.
     *
     * This constructor demonstrates the Composition pattern by accepting
     * Point3D objects. The line is defined by its two endpoints, following
     * the Value Object pattern where the line's identity is its state.
     *
     * The constructor validates inputs and creates an immutable line segment.
     * All subsequent operations will return new Line3D instances rather than
     * modifying this one.
     *
     * @param start the starting point of the line segment
     * @param end the ending point of the line segment
     * @throws IllegalArgumentException if either point is null or if points are identical
     *
     * Example usage:
     * <pre>
     * Point3D p1 = new Point3D(0, 0, 0);
     * Point3D p2 = new Point3D(1, 1, 1);
     * Line3D line = new Line3D(p1, p2);
     * </pre>
     */
    public Line3D(Point3D start, Point3D end) {
        if (start == null || end == null) {
            logger.log(Level.SEVERE, "Attempted to create line with null point(s)");
            throw new IllegalArgumentException("Start and end points cannot be null");
        }

        if (start.equals(end)) {
            logger.log(Level.SEVERE,
                    "Attempted to create line with identical start and end points: {0}",
                    start);
            throw new IllegalArgumentException(
                    "Start and end points must be different (degenerate line not allowed)");
        }

        this.start = start;
        this.end = end;

        logger.log(Level.INFO, "Created new Line3D from {0} to {1}",
                new Object[]{start, end});
    }

    /**
     * Returns the starting point of this line segment.
     *
     * This accessor method demonstrates encapsulation. While the Point3D
     * is immutable, we return the actual reference (not a copy) because
     * Point3D itself cannot be modified.
     *
     * Time Complexity: O(1)
     *
     * @return the start point
     */
    public Point3D getStart() {
        return start;
    }

    /**
     * Returns the ending point of this line segment.
     *
     * Time Complexity: O(1)
     *
     * @return the end point
     */
    public Point3D getEnd() {
        return end;
    }

    /**
     * Calculates the length of this line segment.
     *
     * The length is computed as the Euclidean distance between the start
     * and end points using the distance formula:
     * length = sqrt((x2-x1)² + (y2-y1)² + (z2-z1)²)
     *
     * This method demonstrates composition by delegating to Point3D's
     * distanceTo method rather than reimplementing the distance calculation.
     * This follows the DRY (Don't Repeat Yourself) principle.
     *
     * Algorithm Analysis:
     * - Time Complexity: O(1) - constant time calculation
     * - Space Complexity: O(1) - only local variables used
     *
     * Applications:
     * - Perimeter calculations
     * - Physics (displacement magnitude)
     * - Path planning (distance metrics)
     * - Computer graphics (edge lengths)
     *
     * @return the length of the line segment
     *
     * Example usage:
     * <pre>
     * Line3D line = new Line3D(new Point3D(0, 0, 0), new Point3D(3, 4, 0));
     * double length = line.length(); // Returns 5.0
     * </pre>
     */
    public double length() {
        double len = start.distanceTo(end);
        logger.log(Level.INFO, "Calculated length of line {0}: {1}",
                new Object[]{this, len});
        return len;
    }

    /**
     * Calculates the direction vector of this line (from start to end).
     *
     * The direction vector is computed as: end - start
     * This vector points from the start point toward the end point and
     * has magnitude equal to the line's length.
     *
     * This is a fundamental operation in computational geometry, used for:
     * - Parametric line equations: P(t) = start + t * direction
     * - Parallel/perpendicular line detection
     * - Projection calculations
     * - Intersection algorithms
     *
     * Note: This returns the un-normalized direction vector. For a unit
     * direction vector, use getUnitDirection().
     *
     * Time Complexity: O(1)
     *
     * @return a Point3D representing the direction vector
     *
     * Example usage:
     * <pre>
     * Line3D line = new Line3D(new Point3D(1, 2, 3), new Point3D(4, 6, 8));
     * Point3D direction = line.getDirection(); // Returns Point3D(3, 4, 5)
     * </pre>
     */
    public Point3D getDirection() {
        Point3D direction = end.subtract(start);
        logger.log(Level.INFO, "Calculated direction vector: {0}", direction);
        return direction;
    }

    /**
     * Calculates the unit direction vector of this line.
     *
     * The unit direction vector has magnitude 1 and points in the same
     * direction as the line. This is computed by normalizing the direction
     * vector: unitDirection = direction / |direction|
     *
     * Unit direction vectors are essential in:
     * - Ray tracing (ray directions)
     * - Physics (velocity directions)
     * - Computer graphics (surface normals, light directions)
     * - Geometric calculations requiring direction without magnitude
     *
     * Time Complexity: O(1)
     *
     * @return a Point3D representing the normalized direction vector
     */
    public Point3D getUnitDirection() {
        Point3D unitDir = getDirection().normalize();
        logger.log(Level.INFO, "Calculated unit direction vector: {0}", unitDir);
        return unitDir;
    }

    /**
     * Calculates the midpoint of this line segment.
     *
     * The midpoint is computed as: (start + end) / 2
     * This is the point that divides the line segment into two equal parts.
     *
     * Mathematically, this is the point at parameter t = 0.5 in the
     * parametric equation: P(t) = start + t * (end - start)
     *
     * Applications:
     * - Bisection algorithms
     * - Median calculations
     * - Subdivision algorithms (divide and conquer)
     * - Computer graphics (LOD systems)
     *
     * Time Complexity: O(1)
     *
     * @return the midpoint of the line segment
     *
     * Example usage:
     * <pre>
     * Line3D line = new Line3D(new Point3D(0, 0, 0), new Point3D(10, 10, 10));
     * Point3D mid = line.getMidpoint(); // Returns Point3D(5, 5, 5)
     * </pre>
     */
    public Point3D getMidpoint() {
        Point3D midpoint = start.add(end).multiply(0.5);
        logger.log(Level.INFO, "Calculated midpoint: {0}", midpoint);
        return midpoint;
    }

    /**
     * Gets a point along the line at the specified parameter value.
     *
     * This implements the parametric equation of a line:
     * P(t) = start + t * (end - start)
     *
     * Parameter values:
     * - t = 0: returns the start point
     * - t = 1: returns the end point
     * - 0 < t < 1: returns a point on the line segment
     * - t < 0 or t > 1: returns a point on the infinite line beyond the segment
     *
     * This is fundamental to:
     * - Linear interpolation (LERP)
     * - Animation (tweening)
     * - Sampling points along a path
     * - Subdivision algorithms
     *
     * Time Complexity: O(1)
     *
     * @param t the parameter value (0 = start, 1 = end)
     * @return the point at parameter t
     *
     * Example usage:
     * <pre>
     * Line3D line = new Line3D(new Point3D(0, 0, 0), new Point3D(10, 0, 0));
     * Point3D quarter = line.getPointAt(0.25); // Returns Point3D(2.5, 0, 0)
     * </pre>
     */
    public Point3D getPointAt(double t) {
        Point3D direction = getDirection();
        Point3D point = start.add(direction.multiply(t));

        logger.log(Level.INFO, "Calculated point at t={0}: {1}",
                new Object[]{t, point});

        if (t < 0 || t > 1) {
            logger.log(Level.WARNING,
                    "Parameter t={0} is outside [0,1] range; returning point beyond segment",
                    t);
        }

        return point;
    }

    /**
     * Calculates the shortest distance from a point to this line segment.
     *
     * This is a fundamental computational geometry algorithm with three cases:
     *
     * 1. If the projection of the point onto the infinite line falls on the
     *    segment (0 ≤ t ≤ 1), the distance is the perpendicular distance.
     * 2. If t < 0, the closest point is the start point.
     * 3. If t > 1, the closest point is the end point.
     *
     * Algorithm:
     * 1. Project point onto line: t = ((point - start) · direction) / |direction|²
     * 2. Clamp t to [0, 1] for segment
     * 3. Find closest point: closestPoint = start + t * direction
     * 4. Return distance from point to closestPoint
     *
     * Applications:
     * - Collision detection (point-line proximity)
     * - Robotics (obstacle avoidance)
     * - Computer graphics (selection, snapping)
     * - GIS (nearest road segment)
     *
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     *
     * @param point the point to measure distance from
     * @return the shortest distance from the point to this line segment
     * @throws IllegalArgumentException if point is null
     *
     * Example usage:
     * <pre>
     * Line3D line = new Line3D(new Point3D(0, 0, 0), new Point3D(10, 0, 0));
     * Point3D point = new Point3D(5, 3, 0);
     * double distance = line.distanceToPoint(point); // Returns 3.0
     * </pre>
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Point cannot be null");
        }

        Point3D direction = getDirection();
        Point3D startToPoint = point.subtract(start);

        // Calculate projection parameter
        double directionLengthSquared = direction.dotProduct(direction);
        double t = startToPoint.dotProduct(direction) / directionLengthSquared;

        // Clamp t to [0, 1] for line segment
        t = Math.max(0, Math.min(1, t));

        // Find closest point on segment
        Point3D closestPoint = start.add(direction.multiply(t));

        // Calculate distance
        double distance = point.distanceTo(closestPoint);

        logger.log(Level.INFO,
                "Distance from point {0} to line {1}: {2} (closest point at t={3})",
                new Object[]{point, this, distance, t});

        return distance;
    }

    /**
     * Finds the closest point on this line segment to the given point.
     *
     * This method uses the same projection algorithm as distanceToPoint()
     * but returns the actual point rather than the distance.
     *
     * The algorithm projects the point onto the infinite line, then clamps
     * the result to the segment bounds.
     *
     * Applications:
     * - Snapping in CAD software
     * - Finding nearest point on a path
     * - Projection operations
     * - Collision response (contact points)
     *
     * Time Complexity: O(1)
     *
     * @param point the reference point
     * @return the closest point on this line segment to the given point
     * @throws IllegalArgumentException if point is null
     */
    public Point3D closestPointTo(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Attempted to find closest point to null");
            throw new IllegalArgumentException("Point cannot be null");
        }

        Point3D direction = getDirection();
        Point3D startToPoint = point.subtract(start);

        double directionLengthSquared = direction.dotProduct(direction);
        double t = startToPoint.dotProduct(direction) / directionLengthSquared;

        // Clamp to segment
        t = Math.max(0, Math.min(1, t));

        Point3D closestPoint = start.add(direction.multiply(t));

        logger.log(Level.INFO, "Closest point on line to {0}: {1}",
                new Object[]{point, closestPoint});

        return closestPoint;
    }

    /**
     * Calculates the shortest distance between this line and another line in 3D space.
     *
     * This is a classic computational geometry problem. In 3D, two lines can be:
     * 1. Parallel (including coincident)
     * 2. Intersecting
     * 3. Skew (most common case - non-parallel and non-intersecting)
     *
     * Algorithm for skew lines:
     * The shortest distance occurs along the mutual perpendicular.
     *
     * Given lines:
     * L1: P1(s) = start1 + s * dir1
     * L2: P2(t) = start2 + t * dir2
     *
     * The distance formula is:
     * distance = |((start2 - start1) · (dir1 × dir2))| / |dir1 × dir2|
     *
     * Special cases:
     * - Parallel lines: Use point-to-line distance
     * - Intersecting lines: Distance is 0
     *
     * This demonstrates advanced geometric algorithms and vector mathematics.
     *
     * Applications:
     * - Robotics (collision detection between moving parts)
     * - Computer graphics (line-line proximity tests)
     * - CAD/CAM (clearance calculations)
     * - Physics simulations (particle trajectories)
     * - Network optimization (shortest connection between routes)
     *
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     *
     * @param other the other line
     * @return the shortest distance between the two line segments
     * @throws IllegalArgumentException if other is null
     *
     * Example usage:
     * <pre>
     * // Two skew lines
     * Line3D line1 = new Line3D(new Point3D(0, 0, 0), new Point3D(1, 0, 0));
     * Line3D line2 = new Line3D(new Point3D(0, 1, 1), new Point3D(1, 1, 1));
     * double distance = line1.shortestDistanceTo(line2); // Returns 1.0
     * </pre>
     */
    public double shortestDistanceTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate distance to null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }

        Point3D dir1 = this.getDirection();
        Point3D dir2 = other.getDirection();
        Point3D startDiff = other.start.subtract(this.start);

        // Calculate cross product of directions
        Point3D crossProduct = dir1.crossProduct(dir2);
        double crossMagnitude = crossProduct.magnitude();

        // Check if lines are parallel (cross product ≈ 0)
        if (crossMagnitude < EPSILON) {
            logger.log(Level.INFO, "Lines are parallel; using point-to-line distance");

            // For parallel lines, use point-to-line distance
            // Check multiple points to handle segment boundaries
            double dist1 = this.distanceToPoint(other.start);
            double dist2 = this.distanceToPoint(other.end);
            double dist3 = other.distanceToPoint(this.start);
            double dist4 = other.distanceToPoint(this.end);

            double minDist = Math.min(Math.min(dist1, dist2), Math.min(dist3, dist4));

            logger.log(Level.INFO,
                    "Shortest distance between parallel lines: {0}",
                    minDist);

            return minDist;
        }

        // For skew lines: distance = |(P2 - P1) · (d1 × d2)| / |d1 × d2|
        double distance = Math.abs(startDiff.dotProduct(crossProduct)) / crossMagnitude;

        // However, this is for infinite lines. For segments, we need to check
        // if the closest points actually lie on the segments.

        // Calculate parameters for closest points on infinite lines
        double a = dir1.dotProduct(dir1);
        double b = dir1.dotProduct(dir2);
        double c = dir2.dotProduct(dir2);
        double d = dir1.dotProduct(startDiff);
        double e = dir2.dotProduct(startDiff);

        double denominator = a * c - b * b;

        double s, t;

        if (Math.abs(denominator) < EPSILON) {
            // Lines are parallel (shouldn't reach here, but safety check)
            s = 0;
            t = d / b;
        } else {
            s = (b * e - c * d) / denominator;
            t = (a * e - b * d) / denominator;
        }

        // Clamp parameters to [0, 1] for segments
        s = Math.max(0, Math.min(1, s));
        t = Math.max(0, Math.min(1, t));

        // Calculate actual closest points on segments
        Point3D closestOnThis = this.getPointAt(s);
        Point3D closestOnOther = other.getPointAt(t);

        double segmentDistance = closestOnThis.distanceTo(closestOnOther);

        logger.log(Level.INFO,
                "Shortest distance between lines: {0} (at s={1}, t={2})",
                new Object[]{segmentDistance, s, t});

        return segmentDistance;
    }

    /**
     * Checks if this line is parallel to another line.
     *
     * Two lines are parallel if their direction vectors are parallel,
     * which occurs when their cross product is approximately zero.
     *
     * Mathematical condition: dir1 × dir2 ≈ 0
     *
     * Note: This includes the case of coincident lines (lines that overlap).
     *
     * Applications:
     * - Geometric constraint solving
     * - CAD software (alignment detection)
     * - Computational geometry algorithms
     *
     * Time Complexity: O(1)
     *
     * @param other the other line to check
     * @return true if the lines are parallel, false otherwise
     * @throws IllegalArgumentException if other is null
     */
    public boolean isParallelTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to check parallelism with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }

        Point3D dir1 = this.getUnitDirection();
        Point3D dir2 = other.getUnitDirection();

        Point3D cross = dir1.crossProduct(dir2);
        boolean isParallel = cross.magnitude() < EPSILON;

        logger.log(Level.INFO, "Lines parallel? {0}", isParallel);

        return isParallel;
    }

    /**
     * Checks if this line is perpendicular to another line.
     *
     * Two lines are perpendicular if their direction vectors are perpendicular,
     * which occurs when their dot product is approximately zero.
     *
     * Mathematical condition: dir1 · dir2 ≈ 0
     *
     * Applications:
     * - Geometric constraint solving
     * - Orthogonal projection detection
     * - CAD software (perpendicularity constraints)
     *
     * Time Complexity: O(1)
     *
     * @param other the other line to check
     * @return true if the lines are perpendicular, false otherwise
     * @throws IllegalArgumentException if other is null
     */
    public boolean isPerpendicularTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to check perpendicularity with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }

        Point3D dir1 = this.getDirection();
        Point3D dir2 = other.getDirection();

        double dotProduct = dir1.dotProduct(dir2);
        boolean isPerpendicular = Math.abs(dotProduct) < EPSILON;

        logger.log(Level.INFO, "Lines perpendicular? {0}", isPerpendicular);

        return isPerpendicular;
    }

    /**
     * Calculates the angle between this line and another line in radians.
     *
     * The angle is calculated using the dot product formula:
     * cos(θ) = (dir1 · dir2) / (|dir1| * |dir2|)
     *
     * Returns the acute angle (0 to π/2 radians, or 0 to 90 degrees).
     *
     * Applications:
     * - Mechanical engineering (joint angles)
     * - Computer graphics (angular measurements)
     * - Physics (trajectory analysis)
     *
     * Time Complexity: O(1)
     *
     * @param other the other line
     * @return the angle between the lines in radians (0 to π/2)
     * @throws IllegalArgumentException if other is null
     */
    public double angleBetween(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate angle with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }

        Point3D dir1 = this.getUnitDirection();
        Point3D dir2 = other.getUnitDirection();

        double angle = dir1.angleBetween(dir2);

        // Return the acute angle
        if (angle > Math.PI / 2) {
            angle = Math.PI - angle;
        }

        logger.log(Level.INFO, "Angle between lines: {0} radians", angle);

        return angle;
    }

    /**
     * Reverses the direction of this line (swaps start and end points).
     *
     * Returns a new Line3D with the start and end points swapped.
     * This represents the same geometric line segment but with opposite direction.
     *
     * Applications:
     * - Graph algorithms (bidirectional edges)
     * - Path reversal
     * - Geometric transformations
     *
     * Time Complexity: O(1)
     *
     * @return a new Line3D with reversed direction
     */
    public Line3D reverse() {
        Line3D reversed = new Line3D(end, start);
        logger.log(Level.INFO, "Reversed line direction: {0}", reversed);
        return reversed;
    }

    /**
     * Extends this line segment by the specified factor.
     *
     * The extension occurs symmetrically from the midpoint, maintaining
     * the line's direction and center position.
     *
     * Factor meanings:
     * - factor = 1.0: returns identical line
     * - factor > 1.0: extends the line
     * - factor < 1.0: shrinks the line
     * - factor = 0.0: collapses to midpoint (throws exception)
     *
     * Algorithm:
     * 1. Calculate midpoint
     * 2. Calculate half-direction vector
     * 3. Scale half-direction by factor
     * 4. Create new endpoints: midpoint ± scaled half-direction
     *
     * Time Complexity: O(1)
     *
     * @param factor the scaling factor (must be positive)
     * @return a new Line3D scaled by the factor
     * @throws IllegalArgumentException if factor is not positive
     *
     * Example usage:
     * <pre>
     * Line3D line = new Line3D(new Point3D(0, 0, 0), new Point3D(10, 0, 0));
     * Line3D extended = line.extend(2.0); // Now goes from (-5, 0, 0) to (15, 0, 0)
     * </pre>
     */
    public Line3D extend(double factor) {
        if (factor <= 0) {
            logger.log(Level.SEVERE,
                    "Attempted to extend line by non-positive factor: {0}",
                    factor);
            throw new IllegalArgumentException("Extension factor must be positive");
        }

        Point3D midpoint = getMidpoint();
        Point3D halfDirection = getDirection().multiply(0.5);
        Point3D scaledHalfDirection = halfDirection.multiply(factor);

        Point3D newStart = midpoint.subtract(scaledHalfDirection);
        Point3D newEnd = midpoint.add(scaledHalfDirection);

        Line3D extended = new Line3D(newStart, newEnd);

        logger.log(Level.INFO, "Extended line by factor {0}: {1}",
                new Object[]{factor, extended});

        return extended;
    }

    /**
     * Translates this line by the specified displacement vector.
     *
     * Both endpoints are moved by the same vector, preserving the line's
     * direction and length while changing its position.
     *
     * This is a fundamental geometric transformation used in:
     * - Computer graphics (object movement)
     * - Animation (position updates)
     * - Physics simulations (displacement)
     *
     * Time Complexity: O(1)
     *
     * @param displacement the vector to add to both endpoints
     * @return a new Line3D translated by the displacement
     * @throws IllegalArgumentException if displacement is null
     */
    public Line3D translate(Point3D displacement) {
        if (displacement == null) {
            logger.log(Level.SEVERE, "Attempted to translate by null displacement");
            throw new IllegalArgumentException("Displacement cannot be null");
        }

        Point3D newStart = start.add(displacement);
        Point3D newEnd = end.add(displacement);

        Line3D translated = new Line3D(newStart, newEnd);

        logger.log(Level.INFO, "Translated line by {0}: {1}",
                new Object[]{displacement, translated});

        return translated;
    }

    /**
     * Checks if this line contains the specified point (with tolerance).
     *
     * A point is considered on the line if the distance from the point to
     * the line is less than EPSILON.
     *
     * This is useful for:
     * - Collision detection
     * - Point-line membership testing
     * - Geometric validation
     *
     * Time Complexity: O(1)
     *
     * @param point the point to check
     * @return true if the point lies on this line segment, false otherwise
     * @throws IllegalArgumentException if point is null
     */
    public boolean contains(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Attempted to check if line contains null point");
            throw new IllegalArgumentException("Point cannot be null");
        }

        double distance = distanceToPoint(point);
        boolean containsPoint = distance < EPSILON;

        logger.log(Level.INFO, "Line contains point {0}? {1} (distance: {2})",
                new Object[]{point, containsPoint, distance});

        return containsPoint;
    }

    /**
     * Checks if this line is equal to another object.
     *
     * Two lines are equal if their start and end points are equal.
     * Note: Lines with reversed direction are NOT considered equal by this method.
     *
     * This implements the Value Object pattern: equality based on state.
     *
     * Time Complexity: O(1)
     *
     * @param obj the object to compare with
     * @return true if the objects represent the same line, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        Line3D other = (Line3D) obj;

        boolean isEqual = this.start.equals(other.start) &&
                this.end.equals(other.end);

        logger.log(Level.INFO, "Equality check: {0} equals {1}? {2}",
                new Object[]{this, other, isEqual});

        return isEqual;
    }

    /**
     * Generates a hash code for this line.
     *
     * Ensures that equal lines have equal hash codes, maintaining the
     * hash code contract required for hash-based collections.
     *
     * Time Complexity: O(1)
     *
     * @return the hash code value
     */
    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + start.hashCode();
        result = 31 * result + end.hashCode();
        return result;
    }

    /**
     * Returns a string representation of this line.
     *
     * Format: "Line3D[start -> end]"
     *
     * Time Complexity: O(1)
     *
     * @return a string representation of this line
     */
    @Override
    public String toString() {
        return String.format("Line3D[%s -> %s]", start, end);
    }

    /**
     * Example usage demonstrating the Line3D class capabilities.
     *
     * This main method serves as both documentation and a test suite,
     * showing various operations and geometric calculations.
     */
    public static void main(String[] args) {
        logger.log(Level.INFO, "Starting Line3D demonstration");

        // Create points
        Point3D p1 = new Point3D(0, 0, 0);
        Point3D p2 = new Point3D(10, 0, 0);
        Point3D p3 = new Point3D(0, 5, 0);
        Point3D p4 = new Point3D(0, 15, 0);

        // Create lines
        Line3D line1 = new Line3D(p1, p2);
        Line3D line2 = new Line3D(p3, p4);

        // Basic properties
        System.out.println("Line 1: " + line1);
        System.out.println("Length: " + line1.length());
        System.out.println("Midpoint: " + line1.getMidpoint());
        System.out.println("Direction: " + line1.getDirection());
        System.out.println("Unit Direction: " + line1.getUnitDirection());

        // Point operations
        Point3D testPoint = new Point3D(5, 3, 0);
        System.out.println("\nDistance to point: " + line1.distanceToPoint(testPoint));
        System.out.println("Closest point: " + line1.closestPointTo(testPoint));
        System.out.println("Point at t=0.5: " + line1.getPointAt(0.5));

        // Line-line operations
        System.out.println("\nShortest distance between lines: " +
                line1.shortestDistanceTo(line2));
        System.out.println("Lines parallel? " + line1.isParallelTo(line2));
        System.out.println("Lines perpendicular? " + line1.isPerpendicularTo(line2));
        System.out.println("Angle between lines: " +
                Math.toDegrees(line1.angleBetween(line2)) + " degrees");

        // Transformations
        System.out.println("\nReversed: " + line1.reverse());
        System.out.println("Extended 2x: " + line1.extend(2.0));
        System.out.println("Translated: " +
                line1.translate(new Point3D(5, 5, 5)));

        logger.log(Level.INFO, "Line3D demonstration completed");
    }
}