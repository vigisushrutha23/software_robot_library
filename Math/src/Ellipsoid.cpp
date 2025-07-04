/**
 * @file    Ellipoid.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    May 2025
 * @version 1.0
 * @brief   Source files for the Ellipsoid class.
 *
 * @details An ellipsoid is defined by its centerpoint c, and matrix A such that, for any point p on
 *          it surface it satisfies (p - c)^T A^{-1} (p - c) = 1.
 *
 * @copyright Copyright (c) 2025 Jon Woolfrey
 * 
 * @license GNU General Public License V3
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */
 
#include <Math/Ellipsoid.h>

namespace RobotLibrary { namespace Math {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                         Constructor                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Ellipsoid<Dim>::Ellipsoid(const Eigen::Vector<double,Dim>     &center,
                          const Eigen::Matrix<double,Dim,Dim> &shapeMatrix)
: _center(center),
  _shapeMatrix(shapeMatrix)
{
    _LLT = shapeMatrix.llt();
    
    if (_LLT.info() != Eigen::Success)
    {
        throw std::runtime_error("[ERROR] [ELLIPSOID] Constructor: "
                                 "Shape matrix is not positive definite; Cholesky decomposition failed.");
    }
    _rotationMatrix = Eigen::MatrixXd::Identity(shapeMatrix.rows(),shapeMatrix.rows());
    for(int i = 0; i < shapeMatrix.rows(); i++)
        _axesLengths(i) = sqrt(1/shapeMatrix(i,i));
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                         Constructor                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Ellipsoid<Dim>::Ellipsoid(const Eigen::Vector<double,Dim>     &center,
                          const Eigen::Matrix<double,Dim,Dim> &rotationMatrix, const Eigen::Vector<double,Dim> axesLengths)
: _center(center),
  _axesLengths(axesLengths),
  _rotationMatrix(rotationMatrix)
{
    Eigen::MatrixXd temp_ellipsoid_matrix = Eigen::MatrixXd::Identity(rotationMatrix.rows(),rotationMatrix.rows());
    for(int i = 0; i < rotationMatrix.rows(); i++)
        temp_ellipsoid_matrix(i,i) = 1/pow(axesLengths(i),2);
    _shapeMatrix = _rotationMatrix*temp_ellipsoid_matrix*_rotationMatrix.transpose();

    _LLT = _shapeMatrix.llt();
    
    if (_LLT.info() != Eigen::Success)
    {
        throw std::runtime_error("[ERROR] [ELLIPSOID] Constructor: "
                                 "Shape matrix is not positive definite; Cholesky decomposition failed.");
    }
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                 Compute the distance to a point                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
double
Ellipsoid<Dim>::distance(const Eigen::Vector<double, Dim> &point)
{
    return sqrt(abs(distance_squared(point)));
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                              Compute the squared distance to a point                           //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
double
Ellipsoid<Dim>::distance_squared(const Eigen::Vector<double, Dim> &point)
{
    Eigen::Vector<double, Dim> v = point - _center;
    
    return v.dot(_LLT.solve(v)) - 1.0;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                          Determine if a point is inside the ellipsoid or not                   //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
bool
Ellipsoid<Dim>::is_inside(const Eigen::Vector<double, Dim> &point)
{
    return distance_squared(point) < 0.0;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                  Inverse of the shape matrix multiplied by displacement from center            //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Eigen::Vector<double, Dim>
Ellipsoid<Dim>::inverse_shape_transformed_vector(const Eigen::Vector<double,Dim> &point)
{
    return _LLT.solve(point - _center);
}


  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                              Get inflated ellipsoid  matrix                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Eigen::Matrix<double, Dim, Dim> Ellipsoid<Dim>::get_inflated_ellipsoid_matrix(const double inflation_radius) const
{
    Eigen::MatrixXd temp_ellipsoid_matrix = Eigen::MatrixXd::Identity(_rotationMatrix.rows(),_rotationMatrix.rows());
    for(int i = 0; i < _axesLengths.size(); i++)
        temp_ellipsoid_matrix(i,i) = 1/pow(_axesLengths(i)+inflation_radius,2);
    return _rotationMatrix*temp_ellipsoid_matrix*_rotationMatrix.transpose();
}

} } // namespace
