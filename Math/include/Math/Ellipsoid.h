/**
 * @file    Ellipsoid.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    May 2025
 * @version 1.0
 * @brief   A class that represents an ellipsoid in n-dimensional space.
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

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <Eigen/Dense>                                                                              // Allows Cholesky decomposition

namespace RobotLibrary { namespace Math {

/**
 * @brief A class for representing n-dimensional ellipsoids. Given the center c, and shape matrix A,
 *        an ellipsoid satisfies (p - c)^T * A^-1 * (p - c) = 1 for any point p on its surface.
 */
template <unsigned int Dim>
class Ellipsoid
{
    public:
    
        /**
         * @brief Constructor.
         * @param center The position of the center.
         * @param shapeMatrix A positive-definite matrix.
         */
        Ellipsoid(const Eigen::Vector<double,Dim>     &center,
                  const Eigen::Matrix<double,Dim,Dim> &shapeMatrix);

        /**
         * @brief Constructor.
         * @param center The position of the center.
         * @param rotationMatrix A rotation matrix to orient the ellipse.
         * @param axesLengths Axes lenghts of the ellipsoid
         */
        Ellipsoid(const Eigen::Vector<double,Dim>     &center,
            const Eigen::Matrix<double,Dim,Dim> &rotationMatrix, const Eigen::Vector<double,Dim> axesLengths) ;

                  
        /**
         * @brief Get the distance to a point.
         * @param point What you'd expect.
         * @return What you asked for.
         */
        double
        distance(const Eigen::Vector<double,Dim> &point);
        
        /**
         * @brief Computes the squared distance from the surface of the ellipsoid.
         * @param point What you want the distance to.
         * @return What you asked for.
         */
        double distance_squared(const Eigen::Vector<double,Dim> &point);
        
        /**
         * @brief Used to check if a point is inside the ellipsoid or not.
         * @param point The thing to check.
         * @return Yes or no, true or false.
         */
        bool
        is_inside(const Eigen::Vector<double,Dim> &point);   
        
        /**
         * @brief Get the transformation A^{-1} * (p - c)
         * @param point The point offset from the center.
         * @return The transformed vector.
         */
        Eigen::Vector<double,Dim>
        inverse_shape_transformed_vector(const Eigen::Vector<double,Dim> &point);

        /**
         * @brief Inline function to return centre 'c' of the ellipsoid.
         * @return centre vector 'c' of the ellipsoid.
         */
        inline Eigen::Vector<double,Dim> get_centre() const{ return _center;}

        /**
         * @brief Inline function to return matrix A of the ellipsoid.
         * @return Matrix A of the ellipsoid.
         */
        inline Eigen::Matrix<double, Dim, Dim> get_ellipsod_matrix() const{ return _shapeMatrix;}

        /**
         * @brief Function to calculate and get the elipsoid inflated matrix
         * @param inflation_radius radius to inflate the ellipsoid in m
         */
        Eigen::Matrix<double, Dim, Dim> get_inflated_ellipsoid_matrix(const double inflation_radius) const;   
    
    private:
        
        Eigen::LLT<Eigen::Matrix<double, Dim, Dim>> _LLT;                                           ///< Cholesky decomposition of the shape matrix.

        Eigen::Vector<double, Dim> _center;                                                         ///< Center of the ellipsoid
        
        Eigen::Matrix<double, Dim, Dim> _shapeMatrix;                                               ///< A positive definite matrix describing its shape

        Eigen::Vector<double, Dim> _axesLengths;                                                    ///< Axes lengths of the ellipse

        Eigen::Matrix<double, Dim,Dim> _rotationMatrix;                                                       ///< Rotation matrix to denote orientation of ellipse 

};

} } // namespace

#endif
