use std::ops;

pub type Matrix3 = Matrix<3>;
pub type Matrix4 = Matrix<4>;

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix<const D: usize> {
    m: [[f32; D]; D],
}

impl<const D: usize> Matrix<D> {
    pub fn identity() -> Self {
        let mut m = [[0.0; D]; D];
        for i in 0..D {
            m[i][i] = 1.0;
        }
        Self { m }
    }
    pub fn empty() -> Self {
        Self { m: [[0.0; D]; D] }
    }
}

// Quite epic
impl<const D: usize> ops::Mul<Self> for Matrix<D> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut m = Self::empty();
        for row in 0..D {
            for column in 0..D {
                for i in 0..D {
                    m.m[row][column] += self.m[row][i] * rhs.m[i][column];
                }
            }
        }
        m
    }
}

pub type Vector2 = Vector<2>;
pub type Vector3 = Vector<3>;

#[derive(PartialEq, Debug)]
pub struct Vector<const D: usize> {
    components: [f32; D],
}

impl<const D: usize> Vector<D> {
    pub fn zeroed() -> Self {
        Self {
            components: [0.0; D],
        }
    }
    pub fn magnitude(&self) -> f32 {
        let mut accum: f32 = 0.0;
        for c in self.components {
            accum += c * c;
        }
        accum.sqrt()
    }
    pub fn normalize(mut self) -> Self {
        let mag = self.magnitude();
        for c in &mut self.components {
            *c /= mag;
        }
        self
    }
    pub fn dot(&self, rhs: &Self) -> f32 {
        let mut accum = 0.0;
        for (a, b) in self.components.iter().zip(rhs.components.iter()) {
            accum += a * b;
        }
        accum
    }
}

impl<const D: usize> ops::Mul<f32> for Vector<D> {
    type Output = Self;
    fn mul(mut self, rhs: f32) -> Self::Output {
        for c in &mut self.components {
            *c *= rhs;
        }
        self
    }
}

impl<const D: usize> ops::Div<f32> for Vector<D> {
    type Output = Self;
    fn div(mut self, rhs: f32) -> Self::Output {
        for c in &mut self.components {
            *c /= rhs;
        }
        self
    }
}

impl<const D: usize> ops::Add<Vector<D>> for Vector<D> {
    type Output = Self;
    fn add(mut self, rhs: Vector<D>) -> Self::Output {
        for (a, b) in self.components.iter_mut().zip(rhs.components.iter()) {
            *a += b;
        }
        self
    }
}


impl<const D: usize> ops::AddAssign<Vector<D>> for Vector<D> {
    fn add_assign(&mut self, rhs: Vector<D>) {
        for (a, b) in self.components.iter_mut().zip(rhs.components.iter()) {
            *a += b;
        } 
    }
}

impl<const D: usize> ops::Sub<Vector<D>> for Vector<D> {
    type Output = Self;
    fn sub(mut self, rhs: Vector<D>) -> Self::Output {
        for (a, b) in self.components.iter_mut().zip(rhs.components.iter()) {
            *a -= b;
        }
        self
    }
}

impl<const D: usize> ops::SubAssign<Vector<D>> for Vector<D> {
    fn sub_assign(&mut self, rhs: Vector<D>) {
        for (a, b) in self.components.iter_mut().zip(rhs.components.iter()) {
            *a -= b;
        } 
    }
}

impl Vector<2> {
    #[inline] // I think 'inline' is appropriate here
    pub fn new(x: f32, y: f32) -> Self {
        Self { components: [x, y] }
    }
    #[inline]
    pub fn x(&self) -> f32 {
        self.components[0]
    }
    #[inline]
    pub fn y(&self) -> f32 {
        self.components[1]
    }
}

impl Vector<3> {
    #[inline]
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            components: [x, y, z],
        }
    }
    #[inline]
    pub fn x(&self) -> f32 {
        self.components[0]
    }
    #[inline]
    pub fn y(&self) -> f32 {
        self.components[1]
    }
    #[inline]
    pub fn z(&self) -> f32 {
        self.components[2]
    }
    // TODO: implement 'cross product' (I don't get it :( )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix3_identity_test() {
        let m = Matrix3::identity();
        let expected = Matrix3 {
            m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        };
        assert_eq!(m, expected);
    }

    #[test]
    fn matrix4_identity_test() {
        let m = Matrix4::identity();
        let expected = Matrix4 {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        };
        assert_eq!(m, expected);
    }

    #[test]
    fn matrix_multiplication_test() {
        let m0 = Matrix4 {
            m: [
                [1.0, 2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0, 8.0],
                [9.0, 10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0, 16.0],
            ],
        };
        let m1 = m0.clone();
        let expected = Matrix4 {
            m: [
                [90.0, 100.0, 110.0, 120.0],
                [202.0, 228.0, 254.0, 280.0],
                [314.0, 356.0, 398.0, 440.0],
                [426.0, 484.0, 542.0, 600.0],
            ],
        };
        assert_eq!(m0 * m1, expected);
    }
    #[test]
    fn vector_magnitude_test() {
        let v = Vector2::new(3.0, 4.0);
        assert_eq!(v.magnitude(), 5.0);
    }
    #[test]
    fn vector_normalize_test() {
        let v = Vector2::new(3.0, 4.0);
        let expected = Vector2::new(3.0 / 5.0, 4.0 / 5.0);
        assert_eq!(v.normalize(), expected)
    }
    #[test]
    fn vector_dot_test() {
        let v0 = Vector2::new(1.0, 0.0);
        let v1 = Vector2::new(0.0, 1.0);
        let v2 = Vector2::new(0.5, 0.5);
        assert_eq!(v0.dot(&v0), 1.0);
        assert_eq!(v0.dot(&v1), 0.0);
        assert_eq!(v0.dot(&v2), 0.5);
    }
    #[test]
    fn vector_mulscalar_test() {
        let v = Vector2::new(2.0, 4.0);
        let expected = Vector2::new(4.0, 8.0);
        assert_eq!(v*2.0, expected);
    }
    #[test]
    fn vector_divscalar_test() {
        let v = Vector2::new(2.0, 4.0);
        let expected = Vector2::new(1.0, 2.0);
        assert_eq!(v/2.0, expected);
    }
    #[test]
    fn vector_add_test() {
        let v0 = Vector2::new(1.0, 2.0);
        let v1 = Vector2::new(3.0, 4.0);
        let expected = Vector2::new(4.0, 6.0);
        assert_eq!(v0 + v1, expected);
    }
    #[test]
    fn vector_sub_test() {
        let v0 = Vector2::new(3.0, 4.0);
        let v1 = Vector2::new(1.0, 2.0);
        let expected = Vector2::new(2.0, 2.0);
        assert_eq!(v0 - v1, expected);
    }
    #[test]
    fn vector_addassign_test() {
        let mut v = Vector2::new(2.0, 3.0);
        v += Vector2::new(1.0, 1.0);
        let expected = Vector2::new(3.0, 4.0);
        assert_eq!(v, expected);
    }
    #[test]
    fn vector_subassign_test() {
        let mut v = Vector2::new(2.0, 3.0);
        v -= Vector2::new(1.0, 1.0);
        let expected = Vector2::new(1.0, 2.0);
        assert_eq!(v, expected);
    }
}
