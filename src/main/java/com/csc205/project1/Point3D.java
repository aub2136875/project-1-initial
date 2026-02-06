package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;

    /**
     * Represents a point in three-dimensional Euclidean space.
     *
     * This class demonstrates several fundamental object-oriented design patterns
     * and principles:
     *
     * 1. VALUE OBJECT PATTERN:
     *    - Immutable design: All fields are final and private
     *    - No setters provided - modifications return new instances
     *    - Thread-safe by design due to immutability
     *    - Equality based on state rather than identity
     *
     * 2. ENCAPSULATION:
     *    - Private fields with controlled access through public methods
     *    - Internal state cannot be modified after construction
     *    - Provides abstraction over the coordinate representation
     *
     * 3. SINGLE RESPONSIBILITY PRINCIPLE:
     *    - Focused solely on representing and manipulating 3D points
     *    - Each method has one clear purpose
     *
     * 4. ALGORITHM FOUNDATIONS:
     *    - Euclidean distance calculation (geometric algorithms)
     *    - Vector operations (linear algebra fundamentals)
     *    - Rotation matrices (transformation algorithms)
     *    - Dot and cross products (vector mathematics)
     *
     * 5. DATA STRUCTURE PRINCIPLES:
     *    - Efficient storage using primitive types
     *    - O(1) access to coordinate values
     *    - Immutable structure prevents concurrent modification issues
     *
     * @author Generated for Spring Framework style demonstration
     * @version 1.0
     */
    public class Point3D {

        private static final Logger logger = Logger.getLogger(Point3D.class.getName());

        // Coordinates in 3D space
        private final double x;
        private final double y;
        private final double z;

        // Origin constant for convenience
        public static final Point3D ORIGIN = new Point3D(0, 0, 0);

        /**
         * Constructs a new Point3D with the specified coordinates.
         *
         * This constructor initializes an immutable point in 3D space. Following the
         * Value Object pattern, once created, the point's coordinates cannot be changed.
         * Any transformation operations will return new Point3D instances.
         *
         * @param x the x-coordinate in 3D space
         * @param y the y-coordinate in 3D space
         * @param z the z-coordinate in 3D space
         *
         * Example usage:
         * <pre>
         * Point3D point = new Point3D(3.0, 4.0, 5.0);
         * </pre>
         */
        public Point3D(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
            logger.log(Level.INFO, "Created new Point3D: ({0}, {1}, {2})",
                    new Object[]{x, y, z});
        }

        /**
         * Returns the x-coordinate of this point.
         *
         * This is a simple accessor method demonstrating the encapsulation principle.
         * While the internal representation is hidden, controlled read access is provided.
         *
         * Time Complexity: O(1)
         *
         * @return the x-coordinate value
         */
        public double getX() {
            return x;
        }

        /**
         * Returns the y-coordinate of this point.
         *
         * @return the y-coordinate value
         */
        public double getY() {
            return y;
        }

        /**
         * Returns the z-coordinate of this point.
         *
         * @return the z-coordinate value
         */
        public double getZ() {
            return z;
        }

        /**
         * Calculates the Euclidean distance from this point to another point in 3D space.
         *
         * This method implements the 3D distance formula derived from the Pythagorean theorem:
         * distance = sqrt((x2-x1)² + (y2-y1)² + (z2-z1)²)
         *
         * Algorithm Analysis:
         * - Time Complexity: O(1) - constant time regardless of coordinate values
         * - Space Complexity: O(1) - only local variables used
         *
         * This demonstrates fundamental geometric algorithms used in computer graphics,
         * physics simulations, spatial data structures (like k-d trees), and pathfinding.
         *
         * @param other the target point to calculate distance to
         * @return the Euclidean distance between this point and the other point
         * @throws IllegalArgumentException if other is null
         *
         * Example usage:
         * <pre>
         * Point3D p1 = new Point3D(0, 0, 0);
         * Point3D p2 = new Point3D(3, 4, 0);
         * double distance = p1.distanceTo(p2); // Returns 5.0
         * </pre>
         */
        public double distanceTo(Point3D other) {
            if (other == null) {
                logger.log(Level.SEVERE, "Attempted to calculate distance to null point");
                throw new IllegalArgumentException("Other point cannot be null");
            }

            double dx = this.x - other.x;
            double dy = this.y - other.y;
            double dz = this.z - other.z;

            double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

            logger.log(Level.INFO, "Calculated distance from {0} to {1}: {2}",
                    new Object[]{this, other, distance});

            return distance;
        }

        /**
         * Calculates the Manhattan distance (taxicab distance) to another point.
         *
         * Manhattan distance is the sum of absolute differences of coordinates:
         * distance = |x2-x1| + |y2-y1| + |z2-z1|
         *
         * This metric is useful in grid-based pathfinding algorithms and scenarios
         * where movement is restricted to axis-aligned directions.
         *
         * Time Complexity: O(1)
         *
         * @param other the target point
         * @return the Manhattan distance
         * @throws IllegalArgumentException if other is null
         */
        public double manhattanDistanceTo(Point3D other) {
            if (other == null) {
                logger.log(Level.SEVERE, "Attempted to calculate Manhattan distance to null point");
                throw new IllegalArgumentException("Other point cannot be null");
            }

            double distance = Math.abs(this.x - other.x) +
                    Math.abs(this.y - other.y) +
                    Math.abs(this.z - other.z);

            logger.log(Level.INFO, "Calculated Manhattan distance: {0}", distance);

            return distance;
        }

        /**
         * Calculates the magnitude (length) of the vector from the origin to this point.
         *
         * This is equivalent to the distance from the origin (0,0,0) to this point.
         * In vector terminology, this is the L2 norm or Euclidean norm.
         *
         * Formula: magnitude = sqrt(x² + y² + z²)
         *
         * This is fundamental to vector normalization and many physics calculations.
         *
         * Time Complexity: O(1)
         *
         * @return the magnitude of the position vector
         */
        public double magnitude() {
            double mag = Math.sqrt(x * x + y * y + z * z);
            logger.log(Level.INFO, "Calculated magnitude of {0}: {1}",
                    new Object[]{this, mag});
            return mag;
        }

        /**
         * Returns a normalized version of this point (unit vector with same direction).
         *
         * Normalization creates a vector with magnitude 1 pointing in the same direction.
         * This is essential in computer graphics for lighting calculations, physics
         * for direction vectors, and machine learning for feature scaling.
         *
         * Formula: normalized = (x/mag, y/mag, z/mag) where mag = magnitude()
         *
         * Edge Case: If this point is at the origin (magnitude = 0), returns the origin
         * to avoid division by zero.
         *
         * Time Complexity: O(1)
         *
         * @return a new Point3D representing the normalized vector
         */
        public Point3D normalize() {
            double mag = magnitude();

            if (mag == 0) {
                logger.log(Level.WARNING,
                        "Attempted to normalize zero vector; returning origin");
                return ORIGIN;
            }

            Point3D normalized = new Point3D(x / mag, y / mag, z / mag);
            logger.log(Level.INFO, "Normalized {0} to {1}",
                    new Object[]{this, normalized});

            return normalized;
        }

        /**
         * Rotates this point around the X-axis by the specified angle.
         *
         * This method implements a 3D rotation using the rotation matrix for X-axis:
         * <pre>
         * | 1    0         0      |
         * | 0  cos(θ)  -sin(θ)    |
         * | 0  sin(θ)   cos(θ)    |
         * </pre>
         *
         * Rotation matrices are fundamental to:
         * - 3D graphics and game engines
         * - Robotics and kinematics
         * - Computer vision transformations
         * - Animation systems
         *
         * The algorithm demonstrates linear algebra applications in computer science.
         *
         * Time Complexity: O(1)
         * Space Complexity: O(1)
         *
         * @param angleRadians the rotation angle in radians (positive = counterclockwise
         *                     when looking from positive X toward origin)
         * @return a new Point3D representing the rotated point
         *
         * Example usage:
         * <pre>
         * Point3D p = new Point3D(0, 1, 0);
         * Point3D rotated = p.rotateX(Math.PI / 2); // 90 degrees
         * // Result: approximately (0, 0, 1)
         * </pre>
         */
        public Point3D rotateX(double angleRadians) {
            logger.log(Level.INFO, "Rotating point {0} around X-axis by {1} radians",
                    new Object[]{this, angleRadians});

            double cos = Math.cos(angleRadians);
            double sin = Math.sin(angleRadians);

            double newY = y * cos - z * sin;
            double newZ = y * sin + z * cos;

            return new Point3D(x, newY, newZ);
        }

        /**
         * Rotates this point around the Y-axis by the specified angle.
         *
         * Rotation matrix for Y-axis:
         * <pre>
         * | cos(θ)   0  sin(θ) |
         * |   0      1    0    |
         * | -sin(θ)  0  cos(θ) |
         * </pre>
         *
         * Time Complexity: O(1)
         *
         * @param angleRadians the rotation angle in radians (positive = counterclockwise
         *                     when looking from positive Y toward origin)
         * @return a new Point3D representing the rotated point
         */
        public Point3D rotateY(double angleRadians) {
            logger.log(Level.INFO, "Rotating point {0} around Y-axis by {1} radians",
                    new Object[]{this, angleRadians});

            double cos = Math.cos(angleRadians);
            double sin = Math.sin(angleRadians);

            double newX = x * cos + z * sin;
            double newZ = -x * sin + z * cos;

            return new Point3D(newX, y, newZ);
        }

        /**
         * Rotates this point around the Z-axis by the specified angle.
         *
         * Rotation matrix for Z-axis:
         * <pre>
         * | cos(θ)  -sin(θ)  0 |
         * | sin(θ)   cos(θ)  0 |
         * |   0        0     1 |
         * </pre>
         *
         * Time Complexity: O(1)
         *
         * @param angleRadians the rotation angle in radians (positive = counterclockwise
         *                     when looking from positive Z toward origin)
         * @return a new Point3D representing the rotated point
         */
        public Point3D rotateZ(double angleRadians) {
            logger.log(Level.INFO, "Rotating point {0} around Z-axis by {1} radians",
                    new Object[]{this, angleRadians});

            double cos = Math.cos(angleRadians);
            double sin = Math.sin(angleRadians);

            double newX = x * cos - y * sin;
            double newY = x * sin + y * cos;

            return new Point3D(newX, newY, z);
        }

        /**
         * Adds another point to this point (vector addition).
         *
         * Vector addition is fundamental to:
         * - Physics simulations (force combination)
         * - Animation (displacement)
         * - Pathfinding (position updates)
         *
         * This demonstrates the algebraic structure of vector spaces.
         *
         * Time Complexity: O(1)
         *
         * @param other the point to add
         * @return a new Point3D representing the sum
         * @throws IllegalArgumentException if other is null
         */
        public Point3D add(Point3D other) {
            if (other == null) {
                logger.log(Level.SEVERE, "Attempted to add null point");
                throw new IllegalArgumentException("Other point cannot be null");
            }

            Point3D result = new Point3D(this.x + other.x,
                    this.y + other.y,
                    this.z + other.z);

            logger.log(Level.INFO, "Added {0} + {1} = {2}",
                    new Object[]{this, other, result});

            return result;
        }

        /**
         * Subtracts another point from this point (vector subtraction).
         *
         * Vector subtraction is used to:
         * - Find direction vectors between points
         * - Calculate relative positions
         * - Determine displacement
         *
         * Time Complexity: O(1)
         *
         * @param other the point to subtract
         * @return a new Point3D representing the difference
         * @throws IllegalArgumentException if other is null
         */
        public Point3D subtract(Point3D other) {
            if (other == null) {
                logger.log(Level.SEVERE, "Attempted to subtract null point");
                throw new IllegalArgumentException("Other point cannot be null");
            }

            Point3D result = new Point3D(this.x - other.x,
                    this.y - other.y,
                    this.z - other.z);

            logger.log(Level.INFO, "Subtracted {0} - {1} = {2}",
                    new Object[]{this, other, result});

            return result;
        }

        /**
         * Multiplies this point by a scalar value.
         *
         * Scalar multiplication scales the magnitude while preserving direction.
         * Used extensively in physics for scaling forces, velocities, and accelerations.
         *
         * Time Complexity: O(1)
         *
         * @param scalar the value to multiply by
         * @return a new Point3D representing the scaled point
         */
        public Point3D multiply(double scalar) {
            Point3D result = new Point3D(x * scalar, y * scalar, z * scalar);

            logger.log(Level.INFO, "Multiplied {0} by scalar {1} = {2}",
                    new Object[]{this, scalar, result});

            return result;
        }

        /**
         * Calculates the dot product with another point (treating both as vectors).
         *
         * The dot product formula: a · b = ax*bx + ay*by + az*bz
         *
         * The dot product is fundamental to:
         * - Calculating angles between vectors: cos(θ) = (a · b) / (|a| * |b|)
         * - Projection operations
         * - Determining perpendicularity (dot product = 0)
         * - Lighting calculations in graphics (Lambertian reflection)
         *
         * This demonstrates applications of linear algebra in algorithms.
         *
         * Time Complexity: O(1)
         *
         * @param other the other vector
         * @return the dot product value
         * @throws IllegalArgumentException if other is null
         */
        public double dotProduct(Point3D other) {
            if (other == null) {
                logger.log(Level.SEVERE, "Attempted to calculate dot product with null point");
                throw new IllegalArgumentException("Other point cannot be null");
            }

            double result = this.x * other.x + this.y * other.y + this.z * other.z;

            logger.log(Level.INFO, "Dot product of {0} · {1} = {2}",
                    new Object[]{this, other, result});

            return result;
        }

        /**
         * Calculates the cross product with another point (treating both as vectors).
         *
         * The cross product formula:
         * a × b = (ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx)
         *
         * The cross product is essential for:
         * - Finding perpendicular vectors
         * - Calculating surface normals in 3D graphics
         * - Determining rotation axes
         * - Computing torque in physics
         * - Testing whether points are coplanar
         *
         * Properties:
         * - Result is perpendicular to both input vectors
         * - Magnitude equals area of parallelogram formed by vectors
         * - Anti-commutative: a × b = -(b × a)
         *
         * Time Complexity: O(1)
         *
         * @param other the other vector
         * @return a new Point3D representing the cross product vector
         * @throws IllegalArgumentException if other is null
         */
        public Point3D crossProduct(Point3D other) {
            if (other == null) {
                logger.log(Level.SEVERE, "Attempted to calculate cross product with null point");
                throw new IllegalArgumentException("Other point cannot be null");
            }

            double newX = this.y * other.z - this.z * other.y;
            double newY = this.z * other.x - this.x * other.z;
            double newZ = this.x * other.y - this.y * other.x;

            Point3D result = new Point3D(newX, newY, newZ);

            logger.log(Level.INFO, "Cross product of {0} × {1} = {2}",
                    new Object[]{this, other, result});

            return result;
        }

        /**
         * Calculates the angle between this vector and another vector in radians.
         *
         * Uses the dot product formula: cos(θ) = (a · b) / (|a| * |b|)
         * Therefore: θ = arccos((a · b) / (|a| * |b|))
         *
         * This is used in collision detection, AI behavior (field of view),
         * and determining relative orientations.
         *
         * Edge Case: If either vector has zero magnitude, logs a warning and returns 0.
         *
         * Time Complexity: O(1)
         *
         * @param other the other vector
         * @return the angle in radians (0 to π)
         * @throws IllegalArgumentException if other is null
         */
        public double angleBetween(Point3D other) {
            if (other == null) {
                logger.log(Level.SEVERE, "Attempted to calculate angle with null point");
                throw new IllegalArgumentException("Other point cannot be null");
            }

            double mag1 = this.magnitude();
            double mag2 = other.magnitude();

            if (mag1 == 0 || mag2 == 0) {
                logger.log(Level.WARNING,
                        "Cannot calculate angle with zero-magnitude vector; returning 0");
                return 0;
            }

            double cosAngle = this.dotProduct(other) / (mag1 * mag2);

            // Clamp to [-1, 1] to handle floating point errors
            cosAngle = Math.max(-1.0, Math.min(1.0, cosAngle));

            double angle = Math.acos(cosAngle);

            logger.log(Level.INFO, "Angle between {0} and {1}: {2} radians",
                    new Object[]{this, other, angle});

            return angle;
        }

        /**
         * Checks if this point is equal to another object.
         *
         * Demonstrates the Value Object pattern: equality is based on state
         * (coordinate values) rather than object identity.
         *
         * This override is essential for:
         * - Using Point3D in HashMaps and HashSets
         * - Proper comparison in collections
         * - Unit testing
         *
         * Note: Uses epsilon comparison for floating-point values to handle
         * precision issues.
         *
         * Time Complexity: O(1)
         *
         * @param obj the object to compare with
         * @return true if the objects are equal, false otherwise
         */
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;

            Point3D other = (Point3D) obj;

            final double EPSILON = 1e-10;
            boolean isEqual = Math.abs(this.x - other.x) < EPSILON &&
                    Math.abs(this.y - other.y) < EPSILON &&
                    Math.abs(this.z - other.z) < EPSILON;

            logger.log(Level.INFO, "Equality check: {0} equals {1}? {2}",
                    new Object[]{this, other, isEqual});

            return isEqual;
        }

        /**
         * Generates a hash code for this point.
         *
         * This method must be overridden when equals() is overridden to maintain
         * the hash code contract: equal objects must have equal hash codes.
         *
         * This ensures proper behavior in hash-based collections like HashMap,
         * HashSet, and Hashtable.
         *
         * Time Complexity: O(1)
         *
         * @return the hash code value
         */
        @Override
        public int hashCode() {
            long bits = 17;
            bits = 31 * bits + Double.doubleToLongBits(x);
            bits = 31 * bits + Double.doubleToLongBits(y);
            bits = 31 * bits + Double.doubleToLongBits(z);
            return (int) (bits ^ (bits >>> 32));
        }

        /**
         * Returns a string representation of this point.
         *
         * Provides a human-readable format for debugging and logging.
         * Essential for development and troubleshooting.
         *
         * Format: "Point3D(x, y, z)"
         *
         * Time Complexity: O(1)
         *
         * @return a string representation of this point
         */
        @Override
        public String toString() {
            return String.format("Point3D(%.2f, %.2f, %.2f)", x, y, z);
        }

        /**
         * Example usage demonstrating the Point3D class capabilities.
         *
         * This main method serves as both documentation and a simple test suite,
         * showing how to use the various methods and demonstrating the design patterns
         * in action.
         */
        public static void main(String[] args) {
            logger.log(Level.INFO, "Starting Point3D demonstration");

            // Create points
            Point3D p1 = new Point3D(1, 2, 3);
            Point3D p2 = new Point3D(4, 5, 6);

            // Distance calculations
            System.out.println("Distance: " + p1.distanceTo(p2));
            System.out.println("Manhattan Distance: " + p1.manhattanDistanceTo(p2));

            // Vector operations
            System.out.println("Addition: " + p1.add(p2));
            System.out.println("Dot Product: " + p1.dotProduct(p2));
            System.out.println("Cross Product: " + p1.crossProduct(p2));

            // Rotation
            Point3D rotated = p1.rotateZ(Math.PI / 2);
            System.out.println("Rotated 90° around Z-axis: " + rotated);

            // Normalization
            System.out.println("Normalized: " + p1.normalize());
            System.out.println("Magnitude: " + p1.magnitude());

            logger.log(Level.INFO, "Point3D demonstration completed");
        }
    }
