/**
 * @file    Ellipsoid.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   A class that represents an ellipsoid in n-dimensional space.
 *
 * @details An ellipsoid is defined by its centerpoint c, and matrix A such that, for any point p on
 *          it surface it satisfies (p - c)^T A^{-1} (p - c) = 1. The center is always assumed to be
 *          at zero.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <Math/Shape.h>                                                                             // Base class

#include <Eigen/Dense>                                                                              // Allows Cholesky decomposition

namespace RobotLibrary { namespace Math {

/**
 * @brief A class for representing n-dimensional ellipsoids. Given the center c, and shape matrix A,
 *        an ellipsoid satisfies (p - c)^T * A^-1 * (p - c) = 1 for any point p on its surface.
 * @Note  The center is assumed to be c = 0.
 */
template <unsigned int Dim>
class Ellipsoid : public Shape<Dim>
{
    public:
    
        /**
         * @brief Constructor using a positive-definite matrix.
         * @param center The position of the center.
         * @param shapeMatrix A positive-definite matrix.
         */
        Ellipsoid(const Eigen::Vector<double,Dim>     &center,
                  const Eigen::Matrix<double,Dim,Dim> &shapeMatrix)
        {
            _center = center;
            _shapeMatrix = shapeMatrix;
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
                    
        /**
         * @brief Constructor.
         * @param center The position of the center.
         * @param rotationMatrix A rotation matrix to orient the ellipse.
         * @param axesLengths Axes lenghts of the ellipsoid
         */
        Ellipsoid(const Eigen::Vector<double,Dim>     &center,
        const Eigen::Matrix<double,Dim,Dim> &rotationMatrix, const Eigen::Vector<double,Dim> axesLengths)
        {
            _center = center;
            _rotationMatrix= rotationMatrix;
            _axesLengths = axesLengths;
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

                  
        /**
         * @brief Get a point on the circumference on the ray to the center.
         * @note This method overrides the one in the base class.
         * @param referencePoint An external reference point used for computation.
         * @return What you asked for.
         */
        Eigen::Vector<double,Dim>
        inverse_shape_transformed_vector(const Eigen::Vector<double,Dim> &point);

        /**
         * @brief Inline function to return centre 'c' of the ellipsoid.
         * @return centre vector 'c' of the ellipsoid.
         */
        inline Eigen::Vector<double,Dim> get_centre() const{ return _center;}

        /**
         * @brief Inline function to return axes lenghts of the ellipsoid.
         * @return cvector containing axesLengths.
         */
        inline Eigen::Vector<double,Dim> get_axes_lengths() const{ return _axesLengths;}

        /**
         * @brief Inline function to return matrix A of the ellipsoid.
         * @return Matrix A of the ellipsoid.
         */
        inline Eigen::Matrix<double, Dim, Dim> get_ellipsoid_matrix() const{ return _shapeMatrix;}

        /**
         * @brief Function to calculate and get the elipsoid inflated matrix
         * @param inflation_radius radius to inflate the ellipsoid in m
         * @return The inflated ellipsoid matrix
         */
        Eigen::Matrix<double, Dim, Dim> get_inflated_ellipsoid_matrix(const double inflation_radius) const
        {
            Eigen::MatrixXd temp_ellipsoid_matrix = Eigen::MatrixXd::Identity(_rotationMatrix.rows(),_rotationMatrix.rows()); 
            for(int i = 0; i < _axesLengths.size(); i++)
                temp_ellipsoid_matrix(i,i) = 1/pow(_axesLengths(i)+inflation_radius,2);
            return _rotationMatrix*temp_ellipsoid_matrix*_rotationMatrix.transpose();
        };   
    
        Eigen::Vector<double,Dim>
        point_on_surface(const Eigen::Vector<double,Dim> &referencePoint) const override;
    private:
        
        Eigen::LLT<Eigen::Matrix<double, Dim, Dim>> _LLT;                                           ///< Cholesky decomposition of the shape matrix.
        
        Eigen::Matrix<double, Dim, Dim> _shapeMatrix;                                               ///< A positive definite matrix describing its shape

        Eigen::Vector<double, Dim> _axesLengths;                                                    ///< Axes lengths of the ellipse

        Eigen::Matrix<double, Dim,Dim> _rotationMatrix;                                             ///< Rotation matrix to denote orientation of ellipse 

        Eigen::Vector<double, Dim> _center;                                                         ///< Center of the ellipsoid

};

// Convenience aliases
using Ellipsoid2D = Ellipsoid<2>;
using Ellipsoid3D = Ellipsoid<3>;

} } // namespace

//include <Math/Ellipsoid.tpp>

#endif
