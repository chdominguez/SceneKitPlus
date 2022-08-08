//
//  SceneKitExtensions.swift
//  Atomic
//
//  Created by Christian Dominguez on 18/4/22.
//
import SwiftUI
import SceneKit

///Universal Float (iOS) or CGFloat (macOS) depending on the platform
#if os(macOS)
public typealias UFloat = CGFloat
#elseif os(iOS)
public typealias UFloat = Float
#endif

// Depending on the platform, different color frameworks have to be used
#if os(macOS)
public typealias UColor = NSColor
#elseif os(iOS)
public typealias UColor = UIColor
#endif

public extension Color {
    /// SceneKit better handles either UIColor or NSColor, UColor transforms SwiftUI Color to either depending on the platform. uColor is the universal color for both platforms.
    var uColor: UColor {
        return UColor(self)
    }
}

extension SCNVector3: Equatable {
    public static func == (lhs: SCNVector3, rhs: SCNVector3) -> Bool {
        if lhs.x == rhs.x {
            if lhs.y == rhs.y {
                if  lhs.z == rhs.z {
                    return true
                }
            }
        }
        return false
    }
}

/// Make SCNVector3s to be added and substracted
extension SCNVector3: AdditiveArithmetic {
    
    public static func - (lhs: SCNVector3, rhs: SCNVector3) -> SCNVector3 {
        SCNVector3Make(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z)
    }
    
    public static func + (lhs: SCNVector3, rhs: SCNVector3) -> SCNVector3 {
        SCNVector3Make(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z)
    }
    
    public static var zero: SCNVector3 {
        return SCNVector3(0, 0, 0)
    }
}


extension SCNVector3: VectorArithmetic {
    
    /// Scales and modifies the vector
    public mutating func scale(by rhs: Double) {
        dx *= rhs
        dy *= rhs
        dz *= rhs
    }
    /// The magnitude of the vector
    public var magnitudeSquared: Double {
        return sqrt(dx*dx + dy*dy + dz*dz)
    }
}

/// Useful functions for vectors

extension Double {
    /// Rounds the double to decimal places value
    func rounded(toPlaces places: Int) -> Double {
        let d = pow(10.0, Double(places))
        return (self * d).rounded() / d
    }
}

public extension SCNVector3 {
    
    func roundedCopyTo(n: Int) -> SCNVector3 {
        
        let x = self.dx.rounded(toPlaces: n)
        let y = self.dy.rounded(toPlaces: n)
        let z = self.dz.rounded(toPlaces: n)
        
        return SCNVector3(x, y, z)
    }
    
    mutating func roundTo(n: Int) {
        
        dx = dx.rounded(toPlaces: n)
        dy = dy.rounded(toPlaces: n)
        dz = dz.rounded(toPlaces: n)
    }
    
    // Double counterparts of the x,y and z variables
    var dx: Double { get {Double(x)} set {x = UFloat(newValue)} }
    var dy: Double { get {Double(y)} set {y = UFloat(newValue)} }
    var dz: Double { get {Double(z)} set {z = UFloat(newValue)} }
    
    /// Returns the normalized vector
    func normalized() -> SCNVector3 {
        return SCNVector3Make(UFloat(dx / magnitudeSquared), UFloat(dy / magnitudeSquared), UFloat(dz / magnitudeSquared))
    }
    
    /// Returns the dot product between itself and a second vector.
    func dotProduct(_ v2: SCNVector3) -> Double {
        dx * v2.dx + dy * v2.dy + dz * v2.dz
    }
    
    /// Returns the cross product between itself and a second vector.
    func crossProduct(_ v2: SCNVector3) -> SCNVector3 {
        let nx = y*v2.z - z*v2.y
        let ny = z*v2.x - x*v2.z
        let nz = x*v2.y - y*v2.x
        return SCNVector3Make(nx, ny, nz)
    }
    
    /// Returns the scaled vector without modifying the original vector
    func scaled(by rhs: Double) -> SCNVector3 {
        SCNVector3(dx*rhs, dy*rhs, dz*rhs)
    }
    
    func rotated(by angle: Double, withRespectTo n: SCNVector3) -> SCNVector3 {
        
        let firstTerm = self.scaled(by: cos(angle))
        let secondTerm = (n.crossProduct(self)).scaled(by: sin(angle))
        let thridTerm = n.scaled(by: (1 - cos(angle))*(n.dotProduct(self)))
        
        return firstTerm + secondTerm + thridTerm
        
    }
    
    func getSimd() -> SIMD3<Double> {
        return SIMD3(self)
    }
    
    /// Hability to multiply a vector by a scalar factor in the form v * t
    /// - Parameters:
    ///   - lhs: The vector to be scaled
    ///   - rhs: The scaling factor
    /// - Returns: The scaled vector
    static func * (lhs: SCNVector3, rhs: Double) -> SCNVector3 {
        SCNVector3(lhs.dx * rhs, lhs.dy * rhs, lhs.dz * rhs)
    }
    
    /// Hability to divide a vector by a scalar factor in the form v / t
    /// - Parameters:
    ///   - lhs: The vector to be divided
    ///   - rhs: The dividing factor
    /// - Returns: The divided vector
    static func / (lhs: SCNVector3, rhs: Double) -> SCNVector3 {
        SCNVector3(lhs.dx / rhs, lhs.dy / rhs, lhs.dz / rhs)
    }
    
    
    /// Negates the vector and returns the vector with opposite sign
    /// - Parameter rhs: The vector to be negated
    /// - Returns: The negated vector
    static prefix func - (rhs: SCNVector3) -> SCNVector3 {
        SCNVector3(-rhs.x, -rhs.y, -rhs.z)
    }
    
    
    /// Interpolates between the points a and b by the interpolant t. The parameter t is clamped to the range [0, 1]. This is most commonly used to find a point some fraction of the way along a line between two endpoints (e.g. to move an object gradually between those points).
    /// - Parameters:
    ///   - b: Second position to interpolate to
    ///   - t: Value used to interpolate between a and b
    func lerp(_ b: SCNVector3, _ t: Double) -> SCNVector3 {
        return self + (b - self)*t
    }
      
}

public extension CGFloat {
    /// Initializes a CGFloat with a string
    init?(_ string: String) {
        guard let float = Float(string) else {return nil}
        self = CGFloat(Float(float))
    }
}

public extension SCNQuaternion {
    init(axis: SCNVector3, angle: Double) {
        self = SCNQuaternion(axis.x, axis.y, axis.z, UFloat(angle))
    }
}

public extension SCNMatrix4 {
    
    /// Initializes a SCNMatrix4 from an array of 16 Double values
    /// - Parameter fromArray: Values to initialize the 4x4 matrix in 11, 12, 13... 43, 44 order
    init(_ fromArray: [Double]) {
        
        guard fromArray.count == 16 else { fatalError("SCNMatrix4 require 16 values") }
        
        let m11 = fromArray[0].ufloat
        let m12 = fromArray[1].ufloat
        let m13 = fromArray[2].ufloat
        let m14 = fromArray[3].ufloat
        let m21 = fromArray[4].ufloat
        let m22 = fromArray[5].ufloat
        let m23 = fromArray[6].ufloat
        let m24 = fromArray[7].ufloat
        let m31 = fromArray[8].ufloat
        let m32 = fromArray[9].ufloat
        let m33 = fromArray[10].ufloat
        let m34 = fromArray[11].ufloat
        let m41 = fromArray[12].ufloat
        let m42 = fromArray[13].ufloat
        let m43 = fromArray[14].ufloat
        let m44 = fromArray[15].ufloat
        
        self.init(m11: m11, m12: m12, m13: m13, m14: m14, m21: m21, m22: m22, m23: m23, m24: m24, m31: m31, m32: m32, m33: m33, m34: m34, m41: m41, m42: m42, m43: m43, m44: m44)
    }
    
    static func * (lhs: SCNMatrix4, rhs: Double) -> SCNMatrix4 {
        
        var result = lhs
        
        let urhs = rhs.ufloat
        
        result.m11 *= urhs
        result.m12 *= urhs
        result.m13 *= urhs
        result.m14 *= urhs
        
        result.m21 *= urhs
        result.m22 *= urhs
        result.m23 *= urhs
        result.m24 *= urhs
        
        result.m31 *= urhs
        result.m32 *= urhs
        result.m33 *= urhs
        result.m34 *= urhs
        
        result.m41 *= urhs
        result.m42 *= urhs
        result.m43 *= urhs
        result.m44 *= urhs
        
        return result
        
    }
    
    static func * (lhs: SCNMatrix4, rhs: SCNMatrix4) -> SCNMatrix4 {
        var result = SCNMatrix4()
        result.m11 = lhs.m11*rhs.m11 + lhs.m12*rhs.m21 + lhs.m13*rhs.m31 + lhs.m14*rhs.m41
        result.m21 = lhs.m21*rhs.m11 + lhs.m22*rhs.m21 + lhs.m23*rhs.m31 + lhs.m24*rhs.m41
        result.m31 = lhs.m31*rhs.m11 + lhs.m32*rhs.m21 + lhs.m33*rhs.m31 + lhs.m34*rhs.m41
        result.m41 = lhs.m41*rhs.m11 + lhs.m42*rhs.m21 + lhs.m43*rhs.m31 + lhs.m44*rhs.m41
        result.m12 = lhs.m11*rhs.m12 + lhs.m12*rhs.m22 + lhs.m13*rhs.m32 + lhs.m14*rhs.m42
        result.m22 = lhs.m21*rhs.m12 + lhs.m22*rhs.m22 + lhs.m23*rhs.m32 + lhs.m24*rhs.m42
        result.m32 = lhs.m31*rhs.m12 + lhs.m32*rhs.m22 + lhs.m33*rhs.m32 + lhs.m34*rhs.m42
        result.m42 = lhs.m41*rhs.m12 + lhs.m42*rhs.m22 + lhs.m43*rhs.m32 + lhs.m44*rhs.m42
        result.m13 = lhs.m11*rhs.m13 + lhs.m12*rhs.m23 + lhs.m13*rhs.m33 + lhs.m14*rhs.m43
        result.m23 = lhs.m21*rhs.m13 + lhs.m22*rhs.m23 + lhs.m23*rhs.m33 + lhs.m24*rhs.m43
        result.m33 = lhs.m31*rhs.m13 + lhs.m32*rhs.m23 + lhs.m33*rhs.m33 + lhs.m34*rhs.m43
        result.m43 = lhs.m41*rhs.m13 + lhs.m42*rhs.m23 + lhs.m43*rhs.m33 + lhs.m44*rhs.m43
        result.m14 = lhs.m11*rhs.m14 + lhs.m12*rhs.m24 + lhs.m13*rhs.m34 + lhs.m14*rhs.m44
        result.m24 = lhs.m21*rhs.m14 + lhs.m22*rhs.m24 + lhs.m23*rhs.m34 + lhs.m24*rhs.m44
        result.m34 = lhs.m31*rhs.m14 + lhs.m32*rhs.m24 + lhs.m33*rhs.m34 + lhs.m34*rhs.m44
        result.m44 = lhs.m41*rhs.m14 + lhs.m42*rhs.m24 + lhs.m43*rhs.m34 + lhs.m44*rhs.m44
        return result
    }
}

//MARK: SCNGeometry

public extension SCNGeometry {

    class func triangleFrom(_ vector1: SCNVector3, _ vector2: SCNVector3, _ vector3: SCNVector3) -> SCNGeometry {

        let indices: [Int32] = [0, 1, 2]

        let source = SCNGeometrySource(vertices: [vector1, vector2, vector3])

        let element = SCNGeometryElement(indices: indices, primitiveType: .triangles)

        return SCNGeometry(sources: [source], elements: [element])
    }
}

public extension Double {
    var ufloat: UFloat {
        return UFloat(self)
    }
}
