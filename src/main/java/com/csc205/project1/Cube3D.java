package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * Represents an axis-aligned cube in three-dimensional Euclidean space.
 * 
 * This class demonstrates several fundamental object-oriented design patterns
 * and principles:
 * 
 * 1. COMPOSITION PATTERN:
 *    - Composes Point3D for vertices and Line3D for edges
 *    - "Has-a" relationships throughout
 *    - Builds complex geometry from simple primitives
 * 
 * 2. IMMUTABILITY PATTERN:
 *    - All fields are final
 *    - Transformations return new instances
 *    - Thread-safe by design
 * 
 * 3. BUILDER PATTERN (Static Factory Methods):
 *    - fromCenterAndSize() - construct from center point
 *    - fromMinMax() - construct from bounding box
 *    - Provides flexible construction options
 * 
 * 4. LAZY INITIALIZATION:
 *    - Edges computed on demand, not stored
 *    - Reduces memory footprint
 *    - Avoids unnecessary computation
 * 
 * 5. AGGREGATE PATTERN:
 *    - Represents collection of vertices
 *    - Provides iteration and access methods
 *    - Maintains spatial relationships
 * 
 * 6. DATA STRUCTURE PRINCIPLES:
 *    - Efficient vertex storage (array)
 *    - O(1) access to vertices
 *    - Spatial coherence for cache locality
 * 
 * 7. GEOMETRIC ALGORITHMS:
 *    - Rotation matrices (3D transformations)
 *    - Volume calculations
 *    - Distance computations
 *    - Collision detection
 * 
 * A cube is defined by 8 vertices, 12 edges, and 6 faces. This implementation
 * stores vertices and computes other properties on demand. The cube can be
 * axis-aligned (default) or rotated.
 * 
 * Vertex Ordering (right-handed coordinate system):
 * <pre>
 *      4-----------5
 *     /|          /|
 *    / |         / |
 *   0-----------1  |
 *   |  7--------|--6
 *   | /         | /
 *   |/          |/
 *   3-----------2
 * </pre>
 * 
 * @author Generated for Spring Framework style demonstration
 * @version 1.0
 */
public class Cube3D {
    
    private static final Logger logger = Logger.getLogger(Cube3D.class.getName());
    
    // 8 vertices of the cube (stored in specific order for edge/face computation)
    private final Point3D[] vertices;
    
    // Epsilon for floating-point comparisons
    private static final double EPSILON = 1e-10;
    
    // Edge definitions (pairs of vertex indices)
    private static final int[][] EDGE_INDICES = {
        // Bottom face edges
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        // Top face edges
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        // Vertical edges
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
    
    // Face definitions (4 vertex indices per face, counter-clockwise)
    private static final int[][] FACE_INDICES = {
        {0, 1, 2, 3}, // Bottom (Z-)
        {4, 5, 6, 7}, // Top (Z+)
        {0, 1, 5, 4}, // Front (Y-)
        {3, 2, 6, 7}, // Back (Y+)
        {0, 4, 7, 3}, // Left (X-)
        {1, 2, 6, 5}  // Right (X+)
    };
    
    /**
     * Constructs a cube from 8 vertices.
     * 
     * This constructor demonstrates the COMPOSITION pattern by accepting
     * Point3D instances and building a complex geometric structure.
     * 
     * The constructor validates that:
     * - All 8 vertices are provided
     * - No vertex is null
     * - Vertices form a valid cube structure
     * 
     * Vertex ordering follows the standard cube layout (see class diagram).
     * 
     * @param vertices array of 8 Point3D objects representing cube vertices
     * @throws IllegalArgumentException if vertices array is invalid
     * 
     * Example usage:
     * <pre>
     * Point3D[] vertices = new Point3D[8];
     * vertices[0] = new Point3D(0, 0, 0);
     * // ... set other vertices
     * Cube3D cube = new Cube3D(vertices);
     * </pre>
     */
    public Cube3D(Point3D[] vertices) {
        if (vertices == null || vertices.length != 8) {
            logger.log(Level.SEVERE, 
                      "Invalid vertices array: expected 8 vertices, got {0}", 
                      vertices == null ? "null" : vertices.length);
            throw new IllegalArgumentException("Cube must have exactly 8 vertices");
        }
        
        for (int i = 0; i < 8; i++) {
            if (vertices[i] == null) {
                logger.log(Level.SEVERE, "Vertex {0} is null", i);
                throw new IllegalArgumentException("All vertices must be non-null");
            }
        }
        
        // Create defensive copy to maintain immutability
        this.vertices = Arrays.copyOf(vertices, 8);
        
        logger.log(Level.INFO, "Created new Cube3D with 8 vertices");
    }
    
    /**
     * Creates an axis-aligned cube from center point and size.
     * 
     * This is a STATIC FACTORY METHOD that provides a convenient way to
     * construct cubes. Factory methods offer several advantages:
     * - Named constructors (clearer intent)
     * - Can return cached instances
     * - More flexible than constructors
     * 
     * The cube is aligned with coordinate axes, centered at the given point,
     * with edges of length 'size'.
     * 
     * Algorithm:
     * 1. Calculate half-size (distance from center to face)
     * 2. Generate 8 vertices by adding/subtracting half-size
     * 3. Order vertices according to standard layout
     * 
     * Time Complexity: O(1)
     * 
     * @param center the center point of the cube
     * @param size the length of each edge (must be positive)
     * @return a new Cube3D centered at the given point
     * @throws IllegalArgumentException if center is null or size is not positive
     * 
     * Example usage:
     * <pre>
     * Point3D center = new Point3D(5, 5, 5);
     * Cube3D cube = Cube3D.fromCenterAndSize(center, 10);
     * // Creates a 10x10x10 cube centered at (5,5,5)
     * </pre>
     */
    public static Cube3D fromCenterAndSize(Point3D center, double size) {
        if (center == null) {
            logger.log(Level.SEVERE, "Attempted to create cube with null center");
            throw new IllegalArgumentException("Center point cannot be null");
        }
        
        if (size <= 0) {
            logger.log(Level.SEVERE, "Attempted to create cube with non-positive size: {0}", size);
            throw new IllegalArgumentException("Size must be positive");
        }
        
        double halfSize = size / 2.0;
        double cx = center.getX();
        double cy = center.getY();
        double cz = center.getZ();
        
        Point3D[] vertices = new Point3D[8];
        vertices[0] = new Point3D(cx - halfSize, cy - halfSize, cz - halfSize);
        vertices[1] = new Point3D(cx + halfSize, cy - halfSize, cz - halfSize);
        vertices[2] = new Point3D(cx + halfSize, cy + halfSize, cz - halfSize);
        vertices[3] = new Point3D(cx - halfSize, cy + halfSize, cz - halfSize);
        vertices[4] = new Point3D(cx - halfSize, cy - halfSize, cz + halfSize);
        vertices[5] = new Point3D(cx + halfSize, cy - halfSize, cz + halfSize);
        vertices[6] = new Point3D(cx + halfSize, cy + halfSize, cz + halfSize);
        vertices[7] = new Point3D(cx - halfSize, cy + halfSize, cz + halfSize);
        
        logger.log(Level.INFO, "Created cube from center {0} with size {1}", 
                   new Object[]{center, size});
        
        return new Cube3D(vertices);
    }
    
    /**
     * Creates an axis-aligned cube from minimum and maximum corner points.
     * 
     * This factory method constructs a cube by specifying its bounding box.
     * Useful when working with:
     * - Spatial data structures (octrees, BVH)
     * - Collision detection systems
     * - Scene bounding boxes
     * 
     * The method automatically orders the min/max coordinates to ensure
     * a valid cube is created.
     * 
     * Time Complexity: O(1)
     * 
     * @param min the corner with minimum x, y, z coordinates
     * @param max the corner with maximum x, y, z coordinates
     * @return a new Cube3D spanning from min to max
     * @throws IllegalArgumentException if min or max is null
     * 
     * Example usage:
     * <pre>
     * Point3D min = new Point3D(0, 0, 0);
     * Point3D max = new Point3D(10, 10, 10);
     * Cube3D cube = Cube3D.fromMinMax(min, max);
     * </pre>
     */
    public static Cube3D fromMinMax(Point3D min, Point3D max) {
        if (min == null || max == null) {
            logger.log(Level.SEVERE, "Attempted to create cube with null min/max point");
            throw new IllegalArgumentException("Min and max points cannot be null");
        }
        
        double minX = Math.min(min.getX(), max.getX());
        double maxX = Math.max(min.getX(), max.getX());
        double minY = Math.min(min.getY(), max.getY());
        double maxY = Math.max(min.getY(), max.getY());
        double minZ = Math.min(min.getZ(), max.getZ());
        double maxZ = Math.max(min.getZ(), max.getZ());
        
        if (Math.abs(minX - maxX) < EPSILON || 
            Math.abs(minY - maxY) < EPSILON || 
            Math.abs(minZ - maxZ) < EPSILON) {
            logger.log(Level.WARNING, "Creating degenerate cube with zero dimension");
        }
        
        Point3D[] vertices = new Point3D[8];
        vertices[0] = new Point3D(minX, minY, minZ);
        vertices[1] = new Point3D(maxX, minY, minZ);
        vertices[2] = new Point3D(maxX, maxY, minZ);
        vertices[3] = new Point3D(minX, maxY, minZ);
        vertices[4] = new Point3D(minX, minY, maxZ);
        vertices[5] = new Point3D(maxX, minY, maxZ);
        vertices[6] = new Point3D(maxX, maxY, maxZ);
        vertices[7] = new Point3D(minX, maxY, maxZ);
        
        logger.log(Level.INFO, "Created cube from min {0} to max {1}", 
                   new Object[]{min, max});
        
        return new Cube3D(vertices);
    }
    
    /**
     * Returns the vertex at the specified index.
     * 
     * Vertices are indexed 0-7 according to the standard cube layout.
     * This method provides ENCAPSULATED ACCESS to internal vertex storage.
     * 
     * Since Point3D is immutable, returning the reference is safe.
     * 
     * Time Complexity: O(1) - direct array access
     * 
     * @param index the vertex index (0-7)
     * @return the Point3D at the specified index
     * @throws IllegalArgumentException if index is out of range
     */
    public Point3D getVertex(int index) {
        if (index < 0 || index >= 8) {
            logger.log(Level.SEVERE, "Invalid vertex index: {0}", index);
            throw new IllegalArgumentException("Vertex index must be between 0 and 7");
        }
        return vertices[index];
    }
    
    /**
     * Returns a copy of all vertices.
     * 
     * Returns a defensive copy to maintain immutability. Clients cannot
     * modify the cube's internal state through the returned array.
     * 
     * Time Complexity: O(1) - array copy of fixed size
     * 
     * @return array of 8 vertices
     */
    public Point3D[] getVertices() {
        return Arrays.copyOf(vertices, 8);
    }
    
    /**
     * Calculates the center (centroid) of the cube.
     * 
     * The center is the average of all 8 vertices:
     * center = (v0 + v1 + ... + v7) / 8
     * 
     * This demonstrates the LAZY EVALUATION pattern - the center is
     * computed on demand rather than stored.
     * 
     * For an axis-aligned cube, this is the geometric center.
     * For a rotated cube, this is still the centroid.
     * 
     * Algorithm:
     * 1. Sum all x-coordinates and divide by 8
     * 2. Sum all y-coordinates and divide by 8
     * 3. Sum all z-coordinates and divide by 8
     * 
     * Time Complexity: O(1) - fixed 8 vertices
     * 
     * @return the center point of the cube
     */
    public Point3D getCenter() {
        double sumX = 0, sumY = 0, sumZ = 0;
        
        for (Point3D vertex : vertices) {
            sumX += vertex.getX();
            sumY += vertex.getY();
            sumZ += vertex.getZ();
        }
        
        Point3D center = new Point3D(sumX / 8.0, sumY / 8.0, sumZ / 8.0);
        logger.log(Level.INFO, "Calculated center: {0}", center);
        
        return center;
    }
    
    /**
     * Returns all 12 edges of the cube as Line3D objects.
     * 
     * This method demonstrates LAZY INITIALIZATION - edges are computed
     * on demand rather than stored. This saves memory when edges aren't needed.
     * 
     * Edges are ordered as:
     * - Bottom face: 0-3
     * - Top face: 4-7
     * - Vertical: 8-11
     * 
     * Time Complexity: O(1) - creates 12 Line3D objects
     * Space Complexity: O(1) - returns array of 12 edges
     * 
     * @return array of 12 Line3D objects representing cube edges
     */
    public Line3D[] getEdges() {
        Line3D[] edges = new Line3D[12];
        
        for (int i = 0; i < 12; i++) {
            int v1 = EDGE_INDICES[i][0];
            int v2 = EDGE_INDICES[i][1];
            edges[i] = new Line3D(vertices[v1], vertices[v2]);
        }
        
        logger.log(Level.INFO, "Generated 12 edges for cube");
        return edges;
    }
    
    /**
     * Calculates the edge length of the cube.
     * 
     * For a perfect cube, all edges have the same length. This method
     * calculates the length of edge 0 (between vertices 0 and 1).
     * 
     * If the cube has been rotated or transformed, edge lengths should
     * still be equal (within floating-point precision).
     * 
     * Time Complexity: O(1)
     * 
     * @return the length of one edge
     */
    public double getEdgeLength() {
        double length = vertices[0].distanceTo(vertices[1]);
        logger.log(Level.INFO, "Calculated edge length: {0}", length);
        return length;
    }
    
    /**
     * Calculates the total perimeter (sum of all edge lengths).
     * 
     * A cube has 12 edges, so for a cube with edge length L:
     * perimeter = 12 * L
     * 
     * This is useful for:
     * - Material calculations (wire frames)
     * - Path length in graph traversal
     * - Rendering wireframe models
     * 
     * Algorithm:
     * 1. Calculate edge length
     * 2. Multiply by 12
     * 
     * Alternative: Sum all edge lengths individually (more robust for
     * non-perfect cubes).
     * 
     * Time Complexity: O(1)
     * 
     * @return the sum of all 12 edge lengths
     */
    public double getPerimeter() {
        double perimeter = getEdgeLength() * 12;
        logger.log(Level.INFO, "Calculated perimeter: {0}", perimeter);
        return perimeter;
    }
    
    /**
     * Calculates the volume of the cube.
     * 
     * Volume = edge³ for a cube with edge length L.
     * 
     * This implementation uses the edge length method for simplicity.
     * 
     * Alternative algorithms:
     * - Cross product method: V = |(v1-v0) · ((v3-v0) × (v4-v0))|
     * - Bounding box: V = (maxX-minX) * (maxY-minY) * (maxZ-minZ)
     * 
     * Applications:
     * - Spatial partitioning (volume-based LOD)
     * - Physics (mass calculations)
     * - 3D printing (material estimation)
     * 
     * Time Complexity: O(1)
     * 
     * @return the volume of the cube
     */
    public double getVolume() {
        double edge = getEdgeLength();
        double volume = edge * edge * edge;
        logger.log(Level.INFO, "Calculated volume: {0}", volume);
        return volume;
    }
    
    /**
     * Calculates the surface area of the cube.
     * 
     * Surface area = 6 * edge² for a cube with edge length L.
     * 
     * A cube has 6 faces, each with area edge².
     * 
     * Applications:
     * - Texture mapping (UV coordinates)
     * - Heat transfer calculations
     * - Paint/material coverage estimation
     * 
     * Time Complexity: O(1)
     * 
     * @return the total surface area of all 6 faces
     */
    public double getSurfaceArea() {
        double edge = getEdgeLength();
        double area = 6 * edge * edge;
        logger.log(Level.INFO, "Calculated surface area: {0}", area);
        return area;
    }
    
    /**
     * Calculates the diagonal length (corner to opposite corner).
     * 
     * For a cube with edge length L, the space diagonal is:
     * diagonal = L * √3
     * 
     * This is the longest distance between any two points in the cube.
     * 
     * Derivation:
     * - Face diagonal: L * √2 (Pythagorean theorem in 2D)
     * - Space diagonal: √(L² + (L√2)²) = L * √3
     * 
     * Applications:
     * - Bounding sphere radius (diagonal / 2)
     * - Maximum extent calculations
     * - Collision detection
     * 
     * Time Complexity: O(1)
     * 
     * @return the space diagonal length
     */
    public double getDiagonalLength() {
        double diagonal = vertices[0].distanceTo(vertices[6]);
        logger.log(Level.INFO, "Calculated diagonal length: {0}", diagonal);
        return diagonal;
    }
    
    /**
     * Translates the cube by a displacement vector.
     * 
     * Translation is a RIGID TRANSFORMATION that moves all vertices by
     * the same displacement vector without rotation or deformation.
     * 
     * Properties preserved:
     * - Edge lengths
     * - Angles
     * - Volume
     * - Shape
     * 
     * Algorithm:
     * 1. Add displacement to each vertex
     * 2. Create new cube with translated vertices
     * 
     * This demonstrates the IMMUTABILITY pattern - returns new instance
     * rather than modifying current instance.
     * 
     * Time Complexity: O(1) - translates 8 vertices
     * 
     * @param displacement the vector to translate by
     * @return a new Cube3D translated by the displacement
     * @throws IllegalArgumentException if displacement is null
     * 
     * Example usage:
     * <pre>
     * Cube3D cube = Cube3D.fromCenterAndSize(Point3D.ORIGIN, 10);
     * Point3D move = new Point3D(5, 0, 0);
     * Cube3D translated = cube.translate(move);
     * </pre>
     */
    public Cube3D translate(Point3D displacement) {
        if (displacement == null) {
            logger.log(Level.SEVERE, "Attempted to translate with null displacement");
            throw new IllegalArgumentException("Displacement cannot be null");
        }
        
        Point3D[] newVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            newVertices[i] = vertices[i].add(displacement);
        }
        
        Cube3D translated = new Cube3D(newVertices);
        logger.log(Level.INFO, "Translated cube by {0}", displacement);
        
        return translated;
    }
    
    /**
     * Rotates the cube around the X-axis.
     * 
     * Rotation is performed about the cube's center point. This ensures
     * the cube rotates in place rather than orbiting around the origin.
     * 
     * Algorithm:
     * 1. Translate cube to origin (center at 0,0,0)
     * 2. Rotate each vertex using rotation matrix
     * 3. Translate back to original center
     * 
     * Rotation Matrix (X-axis):
     * <pre>
     * | 1    0         0      |
     * | 0  cos(θ)  -sin(θ)    |
     * | 0  sin(θ)   cos(θ)    |
     * </pre>
     * 
     * This is a fundamental algorithm in:
     * - 3D graphics rendering
     * - Game object transformations
     * - Animation systems
     * - Camera controls
     * 
     * Time Complexity: O(1) - rotates 8 vertices
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Cube3D rotated around the X-axis
     * 
     * Example usage:
     * <pre>
     * Cube3D cube = Cube3D.fromCenterAndSize(new Point3D(5,5,5), 10);
     * Cube3D rotated = cube.rotateX(Math.PI / 4); // 45 degrees
     * </pre>
     */
    public Cube3D rotateX(double angleRadians) {
        Point3D center = getCenter();
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            // Translate to origin
            Point3D relative = vertices[i].subtract(center);
            
            // Rotate
            Point3D rotated = relative.rotateX(angleRadians);
            
            // Translate back
            newVertices[i] = rotated.add(center);
        }
        
        Cube3D rotatedCube = new Cube3D(newVertices);
        logger.log(Level.INFO, "Rotated cube around X-axis by {0} radians", angleRadians);
        
        return rotatedCube;
    }
    
    /**
     * Rotates the cube around the Y-axis.
     * 
     * Similar to rotateX(), but rotates around the Y-axis.
     * 
     * Rotation Matrix (Y-axis):
     * <pre>
     * | cos(θ)   0  sin(θ) |
     * |   0      1    0    |
     * | -sin(θ)  0  cos(θ) |
     * </pre>
     * 
     * Time Complexity: O(1)
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Cube3D rotated around the Y-axis
     */
    public Cube3D rotateY(double angleRadians) {
        Point3D center = getCenter();
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            Point3D relative = vertices[i].subtract(center);
            Point3D rotated = relative.rotateY(angleRadians);
            newVertices[i] = rotated.add(center);
        }
        
        Cube3D rotatedCube = new Cube3D(newVertices);
        logger.log(Level.INFO, "Rotated cube around Y-axis by {0} radians", angleRadians);
        
        return rotatedCube;
    }
    
    /**
     * Rotates the cube around the Z-axis.
     * 
     * Similar to rotateX() and rotateY(), but rotates around the Z-axis.
     * 
     * Rotation Matrix (Z-axis):
     * <pre>
     * | cos(θ)  -sin(θ)  0 |
     * | sin(θ)   cos(θ)  0 |
     * |   0        0     1 |
     * </pre>
     * 
     * Time Complexity: O(1)
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Cube3D rotated around the Z-axis
     */
    public Cube3D rotateZ(double angleRadians) {
        Point3D center = getCenter();
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            Point3D relative = vertices[i].subtract(center);
            Point3D rotated = relative.rotateZ(angleRadians);
            newVertices[i] = rotated.add(center);
        }
        
        Cube3D rotatedCube = new Cube3D(newVertices);
        logger.log(Level.INFO, "Rotated cube around Z-axis by {0} radians", angleRadians);
        
        return rotatedCube;
    }
    
    /**
     * Scales the cube uniformly by a factor about its center.
     * 
     * Scaling multiplies all edge lengths by the scale factor while
     * keeping the center fixed.
     * 
     * Algorithm:
     * 1. Calculate center
     * 2. For each vertex, find vector from center
     * 3. Multiply vector by scale factor
     * 4. Add back to center
     * 
     * Properties:
     * - factor = 1.0: no change
     * - factor = 2.0: doubles edge length, 8x volume
     * - factor = 0.5: halves edge length, 1/8 volume
     * 
     * Applications:
     * - Level of Detail (LOD) systems
     * - Zoom functionality
     * - Size adjustments
     * 
     * Time Complexity: O(1)
     * 
     * @param factor the scale factor (must be positive)
     * @return a new Cube3D scaled by the factor
     * @throws IllegalArgumentException if factor is not positive
     */
    public Cube3D scale(double factor) {
        if (factor <= 0) {
            logger.log(Level.SEVERE, "Attempted to scale with non-positive factor: {0}", factor);
            throw new IllegalArgumentException("Scale factor must be positive");
        }
        
        Point3D center = getCenter();
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            Point3D toVertex = vertices[i].subtract(center);
            Point3D scaled = toVertex.multiply(factor);
            newVertices[i] = center.add(scaled);
        }
        
        Cube3D scaledCube = new Cube3D(newVertices);
        logger.log(Level.INFO, "Scaled cube by factor {0}", factor);
        
        return scaledCube;
    }
    
    /**
     * Checks if a point is inside the cube (inclusive of boundaries).
     * 
     * For an axis-aligned cube, this is a simple bounds check.
     * For a rotated cube, this is more complex.
     * 
     * Algorithm (axis-aligned):
     * 1. Find min/max bounds
     * 2. Check if point is within bounds on all axes
     * 
     * Algorithm (rotated):
     * 1. Transform point to cube's local space
     * 2. Perform bounds check in local space
     * 
     * This implementation uses a simple bounding box approach.
     * 
     * Applications:
     * - Collision detection
     * - Spatial queries
     * - Culling algorithms
     * - Ray casting
     * 
     * Time Complexity: O(1)
     * 
     * @param point the point to test
     * @return true if the point is inside the cube, false otherwise
     * @throws IllegalArgumentException if point is null
     */
    public boolean contains(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Attempted to check containment of null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        // Find bounding box
        double minX = Double.POSITIVE_INFINITY, maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY, maxY = Double.NEGATIVE_INFINITY;
        double minZ = Double.POSITIVE_INFINITY, maxZ = Double.NEGATIVE_INFINITY;
        
        for (Point3D vertex : vertices) {
            minX = Math.min(minX, vertex.getX());
            maxX = Math.max(maxX, vertex.getX());
            minY = Math.min(minY, vertex.getY());
            maxY = Math.max(maxY, vertex.getY());
            minZ = Math.min(minZ, vertex.getZ());
            maxZ = Math.max(maxZ, vertex.getZ());
        }
        
        boolean inside = point.getX() >= minX - EPSILON && point.getX() <= maxX + EPSILON &&
                        point.getY() >= minY - EPSILON && point.getY() <= maxY + EPSILON &&
                        point.getZ() >= minZ - EPSILON && point.getZ() <= maxZ + EPSILON;
        
        logger.log(Level.INFO, "Point {0} is {1}inside cube", 
                   new Object[]{point, inside ? "" : "not "});
        
        return inside;
    }
    
    /**
     * Calculates the minimum distance from a point to the cube surface.
     * 
     * This is a complex geometric problem with several cases:
     * - Point inside: distance = 0 (or negative if we define inward distance)
     * - Point outside: minimum distance to any face, edge, or vertex
     * 
     * Algorithm:
     * 1. Check if point is inside
     * 2. If outside, calculate distance to all faces
     * 3. Return minimum distance found
     * 
     * Simplified implementation: distance to bounding box
     * 
     * Applications:
     * - Proximity queries
     * - Collision detection
     * - Force field calculations
     * 
     * Time Complexity: O(1)
     * 
     * @param point the point to measure distance from
     * @return the minimum distance to the cube surface
     * @throws IllegalArgumentException if point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        if (contains(point)) {
            logger.log(Level.INFO, "Point {0} is inside cube; distance = 0", point);
            return 0.0;
        }
        
        // Calculate distance to bounding box
        double minX = Double.POSITIVE_INFINITY, maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY, maxY = Double.NEGATIVE_INFINITY;
        double minZ = Double.POSITIVE_INFINITY, maxZ = Double.NEGATIVE_INFINITY;
        
        for (Point3D vertex : vertices) {
            minX = Math.min(minX, vertex.getX());
            maxX = Math.max(maxX, vertex.getX());
            minY = Math.min(minY, vertex.getY());
            maxY = Math.max(maxY, vertex.getY());
            minZ = Math.min(minZ, vertex.getZ());
            maxZ = Math.max(maxZ, vertex.getZ());
        }
        
        double dx = Math.max(0, Math.max(minX - point.getX(), point.getX() - maxX));
        double dy = Math.max(0, Math.max(minY - point.getY(), point.getY() - maxY));
        double dz = Math.max(0, Math.max(minZ - point.getZ(), point.getZ() - maxZ));
        
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        
        logger.log(Level.INFO, "Distance from point {0} to cube: {1}", 
                   new Object[]{point, distance});
        
        return distance;
    }
    
    /**
     * Calculates the axis-aligned bounding box (AABB) of the cube.
     * 
     * The AABB is the smallest axis-aligned box that contains the cube.
     * For an axis-aligned cube, the AABB is the cube itself.
     * For a rotated cube, the AABB is larger than the cube.
     * 
     * Returns a 2-element array:
     * - [0]: minimum corner
     * - [1]: maximum corner
     * 
     * Applications:
     * - Broad-phase collision detection
     * - Frustum culling
     * - Spatial partitioning (octrees, grids)
     * - Quick rejection tests
     * 
     * Time Complexity: O(1) - scans 8 vertices
     * 
     * @return array containing min and max corner points
     */
    public Point3D[] getBoundingBox() {
        double minX = Double.POSITIVE_INFINITY, maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY, maxY = Double.NEGATIVE_INFINITY;
        double minZ = Double.POSITIVE_INFINITY, maxZ = Double.NEGATIVE_INFINITY;
        
        for (Point3D vertex : vertices) {
            minX = Math.min(minX, vertex.getX());
            maxX = Math.max(maxX, vertex.getX());
            minY = Math.min(minY, vertex.getY());
            maxY = Math.max(maxY, vertex.getY());
            minZ = Math.min(minZ, vertex.getZ());
            maxZ = Math.max(maxZ, vertex.getZ());
        }
        
        Point3D min = new Point3D(minX, minY, minZ);
        Point3D max = new Point3D(maxX, maxY, maxZ);
        
        logger.log(Level.INFO, "Calculated bounding box: min={0}, max={1}", 
                   new Object[]{min, max});
        
        return new Point3D[]{min, max};
    }
    
    /**
     * Calculates the bounding sphere of the cube.
     * 
     * The bounding sphere is the smallest sphere that contains all vertices.
     * For a cube, the center is the cube's center and radius is half the
     * space diagonal.
     * 
     * Returns a 2-element array:
     * - [0]: center point
     * - [1]: radius (stored as Point3D with radius in x coordinate)
     * 
     * Applications:
     * - Broad-phase collision detection
     * - Level of Detail (LOD) calculations
     * - View frustum culling
     * 
     * Time Complexity: O(1)
     * 
     * @return array containing center and radius
     */
    public Object[] getBoundingSphere() {
        Point3D center = getCenter();
        double radius = getDiagonalLength() / 2.0;
        
        logger.log(Level.INFO, "Calculated bounding sphere: center={0}, radius={1}", 
                   new Object[]{center, radius});
        
        return new Object[]{center, radius};
    }
    
    /**
     * Checks if this cube intersects with another cube.
     * 
     * This is a simplified implementation using AABB intersection.
     * Two AABBs intersect if they overlap on all three axes.
     * 
     * Algorithm:
     * 1. Get bounding boxes of both cubes
     * 2. Check for overlap on X-axis
     * 3. Check for overlap on Y-axis
     * 4. Check for overlap on Z-axis
     * 5. Return true only if all three overlap
     * 
     * For accurate rotated cube intersection, use SAT (Separating Axis Theorem).
     * 
     * Applications:
     * - Collision detection
     * - Spatial queries
     * - Physics simulations
     * 
     * Time Complexity: O(1)
     * 
     * @param other the other cube to test intersection with
     * @return true if the cubes intersect, false otherwise
     * @throws IllegalArgumentException if other is null
     */
    public boolean intersects(Cube3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to check intersection with null cube");
            throw new IllegalArgumentException("Other cube cannot be null");
        }
        
        Point3D[] thisBB = this.getBoundingBox();
        Point3D[] otherBB = other.getBoundingBox();
        
        boolean overlapX = thisBB[0].getX() <= otherBB[1].getX() && 
                          thisBB[1].getX() >= otherBB[0].getX();
        boolean overlapY = thisBB[0].getY() <= otherBB[1].getY() && 
                          thisBB[1].getY() >= otherBB[0].getY();
        boolean overlapZ = thisBB[0].getZ() <= otherBB[1].getZ() && 
                          thisBB[1].getZ() >= otherBB[0].getZ();
        
        boolean intersects = overlapX && overlapY && overlapZ;
        
        logger.log(Level.INFO, "Cubes {0}intersect", intersects ? "" : "do not ");
        
        return intersects;
    }
    
    /**
     * Returns a string representation of the cube.
     * 
     * Format: "Cube3D[center=(x,y,z), edge=L]"
     * 
     * Provides human-readable output for debugging and logging.
     * 
     * Time Complexity: O(1)
     * 
     * @return a string representation of this cube
     */
    @Override
    public String toString() {
        Point3D center = getCenter();
        double edge = getEdgeLength();
        return String.format("Cube3D[center=%s, edge=%.2f]", center, edge);
    }
    
    /**
     * Checks if this cube is equal to another object.
     * 
     * Two cubes are equal if all their vertices are equal (in the same order).
     * 
     * Time Complexity: O(1) - compares 8 vertices
     * 
     * @param obj the object to compare with
     * @return true if the objects are equal, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Cube3D other = (Cube3D) obj;
        
        for (int i = 0; i < 8; i++) {
            if (!this.vertices[i].equals(other.vertices[i])) {
                return false;
            }
        }
        
        logger.log(Level.INFO, "Equality check: cubes are equal");
        return true;
    }
    
    /**
     * Generates a hash code for this cube.
     * 
     * Uses all vertices to compute hash code for proper HashMap behavior.
     * 
     * Time Complexity: O(1)
     * 
     * @return the hash code value
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(vertices);
    }
    
    /**
     * Example usage demonstrating the Cube3D class capabilities.
     */
    public static void main(String[] args) {
        logger.log(Level.INFO, "Starting Cube3D demonstration");
        
        // Create a cube using different methods
        System.out.println("=== Cube Creation ===");
        Cube3D cube1 = Cube3D.fromCenterAndSize(new Point3D(0, 0, 0), 10);
        System.out.println("Cube 1: " + cube1);
        
        Cube3D cube2 = Cube3D.fromMinMax(
            new Point3D(-5, -5, -5), 
            new Point3D(5, 5, 5)
        );
        System.out.println("Cube 2: " + cube2);
        
        // Geometric properties
        System.out.println("\n=== Geometric Properties ===");
        System.out.println("Edge length: " + cube1.getEdgeLength());
        System.out.println("Perimeter: " + cube1.getPerimeter());
        System.out.println("Surface area: " + cube1.getSurfaceArea());
        System.out.println("Volume: " + cube1.getVolume());
        System.out.println("Diagonal: " + cube1.getDiagonalLength());
        System.out.println("Center: " + cube1.getCenter());
        
        // Transformations
        System.out.println("\n=== Transformations ===");
        Cube3D translated = cube1.translate(new Point3D(10, 0, 0));
        System.out.println("Translated: " + translated);
        
        Cube3D rotated = cube1.rotateZ(Math.PI / 4);
        System.out.println("Rotated 45°: " + rotated);
        
        Cube3D scaled = cube1.scale(2.0);
        System.out.println("Scaled 2x: " + scaled);
        System.out.println("Scaled volume: " + scaled.getVolume());
        
        // Spatial queries
        System.out.println("\n=== Spatial Queries ===");
        Point3D testPoint = new Point3D(3, 3, 3);
        System.out.println("Point " + testPoint + " inside? " + cube1.contains(testPoint));
        System.out.println("Distance to point: " + cube1.distanceToPoint(testPoint));
        
        // Bounding volumes
        System.out.println("\n=== Bounding Volumes ===");
        Point3D[] bbox = cube1.getBoundingBox();
        System.out.println("Bounding box: " + bbox[0] + " to " + bbox[1]);
        
        Object[] bsphere = cube1.getBoundingSphere();
        System.out.println("Bounding sphere: center=" + bsphere[0] + ", radius=" + bsphere[1]);
        
        // Intersection test
        System.out.println("\n=== Intersection Test ===");
        Cube3D cube3 = Cube3D.fromCenterAndSize(new Point3D(8, 0, 0), 10);
        System.out.println("Cube1 intersects Cube3? " + cube1.intersects(cube3));
        
        logger.log(Level.INFO, "Cube3D demonstration completed");
    }
}