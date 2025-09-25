/**
 * @file    DifferentialDrivePredictive.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    May 2025
 * @version 1.0
 * @brief   Source files for the MPC class.
 *
 * @details This class solves the predictive control problem for differential drive robot over a
 *          finite number of steps. It uses quadratic functions for both the final step, and all
 *          intermediate steps. The weighting on the final & intermediate pose errors are used as
 *          constructor arguments, whereas the weighting on intermediate control values is based on
 *          the robot's mass & inertia in the RobotLibary::Model::DifferentialDrive class.
 * 
 * @copyright Copyright (c) 2025 Jon Woolfrey
 * 
 * @license GNU General Public License V3
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 * @see http://github.com/Woolfrey/software_simple_qp for info on the QP solver.
 */
 
#include <Control/DifferentialDrivePredictive.h>

namespace RobotLibrary { namespace Control {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                            Constructor                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////
DifferentialDrivePredictive::DifferentialDrivePredictive(RobotLibrary::Model::DifferentialDriveParameters &modelParameters,
                                                         RobotLibrary::Control::DifferentialDrivePredictiveParameters &controlParameters,
                                                         SolverOptions<double> &solverOptions)
: DifferentialDriveBase(controlParameters.controlFrequency,
                        controlParameters.minimumSafeDistance,
                        modelParameters),
  QPSolver<double>(solverOptions),
 _numberOfRecursions(controlParameters.numberOfRecursions),
 _predictionSteps(controlParameters.predictionSteps),
 _threshold(controlParameters.maximumControlStepNorm)
{
    // Ensure weighting matrices are positive definite:
    std::string message;
    if (not RobotLibrary::Math::is_positive_definite(controlParameters.poseErrorWeight, message))
    {
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE MPC] Constructor: "
                                    "Initial pose error weight matrix is not positive definite: " + message);
    }
    
    // Set size of vectors
    _predictedStates.resize(_predictionSteps + 1);                                                  // 1 for current state + N for predicted states
    _poseErrorWeight.resize(_predictionSteps);
    _controlWeight.resize(_predictionSteps);

    // --- Pose Error Weights (Normalized Exponential) ---
    double exponent = controlParameters.exponent;

    // Precompute denominator for normalization
    double denominator = 0.0;
    for (int j = 0; j < _predictionSteps; ++j)
    {
        denominator += std::exp(exponent * j);
    }

    // Assign normalized exponential pose weights
    for (int j = 0; j < _predictionSteps; ++j)
    {
        double scalar = std::exp(exponent * j) / denominator;
        _poseErrorWeight[j] = scalar * controlParameters.poseErrorWeight;
    }

    // --- Control Weights (Exponentially Scaled Inertia Matrix) ---
    Eigen::Matrix2d inertiaMatrix;
    inertiaMatrix << _mass,  0.0,
                      0.0, _inertia;

    std::vector<double> expWeights(_predictionSteps);

    // Compute exponential profile
    for (int i = 0; i < _predictionSteps; ++i)
    {
        expWeights[i] = std::exp(exponent * i);
    }

    // Anchor at R_0 = M or R_{N-1} = M depending on direction of growth
    double scale = (exponent >= 0.0)
                 ? 1.0 / expWeights.back()   // Ensure R_{N-1} = M
                 : 1.0 / expWeights.front(); // Ensure R_0 = M

    // Assign scaled control effort weights
    for (int i = 0; i < _predictionSteps; ++i)
    {
        double weight = expWeights[i] * scale;
        _controlWeight[i] = weight * inertiaMatrix;
    }
    
    // Set up the constraint matrices in advance to save time:

    _controlConstraintMatrix << -1.0,  0.0,
                                 0.0, -1.0,
                                 1.0,  0.0,
                                 0.0,  1.0;
                                 
    _obstacleConstraintMatrix.resize(0,2);                                                         
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                        Update the state                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
void
DifferentialDrivePredictive::update_state(const RobotLibrary::Model::Pose2D &pose,
                                          const Eigen::Vector2d &velocity,
                                          const Eigen::Matrix3d &covariance)
{
    RobotLibrary::Model::DifferentialDrive::update_state(pose, velocity, covariance);               // Update the underlying model
    
    // Transfer these from the base class so we can access them via index in the
    // backward / forward recursions
    _predictedStates[0].pose       = _pose;
    _predictedStates[0].velocity   = _velocity;
    _predictedStates[0].covariance = _covariance;
    
    // Shift the predicted states backward
    for (int i = 1; i < _predictionSteps - 1; ++i)
    {
        _predictedStates[i] = _predictedStates[i+1];
    }
} 

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                               Solve the trajectory tracking problem                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector2d
DifferentialDrivePredictive::track_trajectory(const std::vector<RobotLibrary::Model::DifferentialDriveState>   &desiredStates,
                                              const std::vector<std::vector<RobotLibrary::Math::Ellipsoid<2>>> &obstacles)
{
    // Ensure inputs are sound
    if (desiredStates.size() != _predictionSteps + 1)
    {
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                    "This controller requires N + 1 = " + std::to_string(_predictionSteps+1) + " "
                                    "desired states for the trajectory tracking, but received "
                                    + std::to_string(desiredStates.size()) + ".");
    }
    
    // Run the optimisations
    for (int i = 0; i < _numberOfRecursions; ++i)
    {
        double largestStepChange = 0.0;                                                             // Store largest step change in control for this recursion

        Eigen::Vector3d costateVector;                                                              // i.e. Lagrange multipliers
        
        // Backwards recursions
        for (int j = _predictionSteps-1; j >= 0; --j)
        {    
            Eigen::Vector3d poseError = _predictedStates[j+1].pose.error(desiredStates[j+1].pose);  // Error at step j+1 is affected by control input at step j
            if (j == _predictionSteps - 1)
            {
                costateVector = _poseErrorWeight[j] * poseError;                                    // Only need to evaluate costate vector at final steps   
            }
            else
            {
                // Values used in this scope           
                RobotLibrary::Model::Pose2D currentPose = _predictedStates[j].pose;
                Eigen::Vector2d currentVelocity = _predictedStates[j].velocity;
                Eigen::Vector2d desiredVelocity = desiredStates[j].velocity;
                Eigen::Matrix3d K = _poseErrorWeight[j];
                Eigen::Matrix2d M = _controlWeight[j];
                double angle = currentPose.angle();                                                 // Used in multiple places
                RobotLibrary::Model::Pose2D desiredPose = desiredStates[j].pose;
                RobotLibrary::Model::Pose2D desiredNextPose =  RobotLibrary::Model::DifferentialDrive::predicted_pose(desiredStates[j].pose,
                    desiredStates[j].velocity,
                    _controlFrequency);

                double angle_desired = desiredPose.angle();                                                 // Used in multiple places

                Eigen::Vector3d poseError_curr = currentPose.error(desiredPose);
 
                // Partial derivative of state propagation w.r.t. configuration                                                    
                Eigen::Matrix<double,3,3> dfdx = configuration_jacobian(currentPose, currentVelocity, _controlFrequency);
                Eigen::Matrix<double,3,3> dfdx_nom = configuration_jacobian(desiredPose, desiredVelocity, _controlFrequency);

                    
                // Partial derivative of state propagation w.r.t. control.
                Eigen::Matrix<double,3,2> dfdu = control_jacobian(currentPose, _controlFrequency);
                Eigen::Matrix<double,3,2> dfdudx = Eigen::MatrixXd::Zero(3,2);
                dfdudx (0,0) = -sin(angle)/_controlFrequency;
                dfdudx (1,0) = cos(angle)/_controlFrequency;

                //dfdu(2,1) = 1.0; // THIS WORKS BETTER?
                Eigen::Vector2d del_u = desiredVelocity - currentVelocity;
                Eigen::Vector3d del_x = currentPose.error(desiredPose);
                Eigen::Vector3d E =   del_x + dfdu * del_u;

                // Partial derivative of state propagation w.r.t. configuration   
                Eigen::Vector<double,2> dLdu =  -del_u.transpose() * M  - poseError.transpose() * K * dfdu + costateVector.transpose()* dfdu;
                                                                      
                // Second derivative of Lagrangian w.r.t. control (i.e. Hessian)
                Eigen::Matrix<double,2,2> d2Ldu2 = M + dfdu.transpose() * K * dfdu;
                                                 
                d2Ldu2(0,0) += 1e-06;                                                               // Add some damping to ensure stability
                d2Ldu2(1,1) += 1e-06;   
                
                // Mixed derivatives of Lagrangian 
                Eigen::Matrix<double,2,3> d2Ldudx = dfdu.transpose() * K * dfdx ;
                d2Ldudx(0,2) += (costateVector.transpose() * dfdudx.col(0));
                d2Ldudx(0,2) -= (poseError.transpose()* K * dfdudx.col(0) );                                                          
                // Set up the constraint for the control input du
                RobotLibrary::Model::Limits linear, angular;
                
                compute_control_limits(linear, angular, currentVelocity);                           // Limits are computed w.r.t. current velocity
                
                _controlConstraintVector << currentVelocity[0] - linear.lower,                      // -dv <= v - v_min
                                            currentVelocity[1] - angular.lower,                     // -dw <= w - w_min
                                            linear.upper  - currentVelocity[0],                     //  dv <= v_max - v
                                            angular.upper - currentVelocity[1];                     //  dw <= w_max - ws
                // Set up the constraints for the obstacles
                /*
                for (int k = 0; k < obstacles.size(); ++k)
                {
                    if (obstacles[k].size() != _predictionSteps + 1)
                    {
                        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                                    "This controller has N + 1 = " + std::to_string(_predictionSteps+1) + " control steps "
                                                    "but obstacle #" + std::to_string(k+1) + " had " + std::to_string(obstacles[k].size()) + " "
                                                    "predicted positions.");
                    }
                    
                    _obstacleConstraintMatrix.row(k) = ...
                    _obstacleConstraintVector.row(k) = ...
                }
                */
               std::cout<<"\n=====Current Pose & Velocity "<<currentPose.translation()(0)<<"\t"<< currentPose.translation()(1)<<"\t"<<currentVelocity(0)<<"\t"<<currentVelocity(1)<<" ====\n";

               Eigen::Vector3d dx = currentPose.error(desiredStates[j].pose);
               _obstacleConstraintMatrix.resize(obstacles.size()*3,2);
               _obstacleConstraintVector.resize(obstacles.size()*3);
                Eigen::MatrixXd P = Eigen::MatrixXd::Zero(2,3); //Selection matrix
                P(0,0) = P(1,1) = 1.0;
                for (int k = 0; k < obstacles.size(); ++k)
                {
                    if (obstacles[k].size() != _predictionSteps)
                    {
                        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                                    "This controller has N = " + std::to_string(_predictionSteps+1) + " control steps "
                                                    "but obstacle #" + std::to_string(k+1) + " had " + std::to_string(obstacles[k].size()) + " "
                                                    "predicted positions.");
                    }
                    Eigen::Vector2d obs_centre = obstacles[k][j].get_centre();
                    
                    Eigen::Vector3d delc = {0.0, 0.0,0.0}; //Stationary object for now
                    for(int l = 0; l<3; l++)
                    {
                        Eigen::Vector3d v_temp_curr  = {0.0,0.0,0.0}; 
                        v_temp_curr.head(2) = currentPose.translation();
                        v_temp_curr(0) += _robotFootprint(l) * cos(angle-M_PI/2) ;
                        v_temp_curr(1) += _robotFootprint(l) * sin(angle-M_PI/2) ;
                        v_temp_curr(2) = angle;
                       
                        Eigen::Vector3d v_temp_desired  = {0.0,0.0,0.0}; 
                        v_temp_desired.head(2) = currentPose.translation();
                        v_temp_desired(0) += _robotFootprint(l) * cos(angle-M_PI/2) + currentVelocity(0)*cos(angle)/_controlFrequency ;
                        v_temp_desired(1) += _robotFootprint(l) * sin(angle-M_PI/2) + currentVelocity(0)*sin(angle)/_controlFrequency ;
                        v_temp_desired(2) = angle + currentVelocity(1)/_controlFrequency;
                       
                        const double temp_radius = _robotRadii(l);
                       
                        Eigen::Matrix2d ellipsoid_shape = obstacles[k][j].get_inflated_ellipsoid_matrix(temp_radius);
                                              

                        Eigen::Vector2d V_desired = P*v_temp_desired - obs_centre ;  //still assuming stationary object
                        Eigen::Vector2d V_curr = P*v_temp_curr - obs_centre ; 

                        

    
                        //Eigen::Vector3d dhdx = V.transpose() * dvdx;
                        //Eigen::Vector2d dhdc = V_curr; 

                        

                        double m = pow(V_desired.transpose() * ellipsoid_shape * V_desired , 0.5);
                        Eigen::Vector2d b;
                        b << cos(v_temp_desired(2)), sin(v_temp_desired(2));
                        Eigen::Matrix2d dbdu = Eigen::Matrix2d::Zero();
                        dbdu(0,1) = -sin(v_temp_desired(2))/_controlFrequency;
                        dbdu(1,1) = cos(v_temp_desired(2))/_controlFrequency;
                        double m_curr = pow(V_curr.transpose() * ellipsoid_shape * V_curr,0.5);
                        

                        Eigen::Matrix<double,2,3> dbdx = Eigen::MatrixXd::Zero(2,3);
                        dbdx(0,2) = -sin(v_temp_desired(2));
                        dbdx(1,2) = cos(v_temp_desired(2));
                        
                        Eigen::Vector2d r = V_desired - V_desired/m;
                        Eigen::Vector2d r_curr = V_curr - V_curr/m_curr;

                        Eigen::Vector2d b_curr;
                        b_curr << cos(v_temp_curr(2)), sin(v_temp_curr(2));

                        Eigen::Matrix2d drdu = (P * dfdu) - (P * dfdu)/m  + (V_desired *V_desired.transpose()*ellipsoid_shape* (P * dfdu) )/pow(m,3);
                        Eigen::Matrix<double,2,3>  drdx = (P * dfdx) - (P * dfdx)/m  + (V_desired *V_desired.transpose()*ellipsoid_shape)*(P * dfdx)/pow(m,3) ;

                        double h_desired =  V_desired.transpose() * ellipsoid_shape * V_desired + b.dot((r/r.norm())) -1.0;



                        Eigen::Vector2d dhdu = 2 * V_desired.transpose()*ellipsoid_shape*P*dfdu + (1/r.norm())*b.transpose()*((Eigen::Matrix2d::Identity() - r*r.transpose()/pow(r.norm(),2))*drdu) + (r.transpose()/r.norm())*dbdu;
                        Eigen::Vector3d dhdx = 2 * V_desired.transpose()*ellipsoid_shape*P*dfdx + (1/r.norm())*b.transpose()*((Eigen::Matrix2d::Identity() - r*r.transpose()/pow(r.norm(),2))*drdx) + (r.transpose()/r.norm())*dbdx;

                        double h_curr =  V_curr.transpose() * ellipsoid_shape * V_curr + b_curr.dot((r_curr/r_curr.norm())) -1.0;
                        std::cout<<"\n Direction Dot product "<<b_curr.dot((r_curr/r_curr.norm()));

                        double alpha = 5.0;
                        double epsilon = 0.0;
                        _obstacleConstraintMatrix.row(k*3+l) = -dhdu.transpose() ;

                        _obstacleConstraintVector(k*3+l) =  h_desired + epsilon + dhdx.transpose()*dx ;
                        if (l==1)
                            dx = v_temp_desired - v_temp_curr;
                        std::cout<<"\n==== Obstacle constraint "<<_obstacleConstraintVector(k*3+l)<<"\t"<< _obstacleConstraintMatrix.row(k*3+l)<<"\t"<<_obstacleConstraintVector(k*3+l)-_obstacleConstraintMatrix.row(k*3+l)*del_u;
                        std::cout<<"\n==== Dfdu \n" << dfdu;
                        std::cout<<"\n==== H_CURR =  "<<  h_desired;
                    }
                }
                        
                // Combine the constraints
                //unsigned int numRows = _controlConstraintVector.size() + _obstacleConstraintVector.size();
                unsigned int numRows =  _obstacleConstraintVector.size();
                _constraintMatrix.resize(numRows, 2);
                //_constraintMatrix.block(0,0,4,2)         = _controlConstraintMatrix;
                _constraintMatrix.block(0,0,numRows,2) = _obstacleConstraintMatrix;
                
                _constraintVector.resize(numRows);
                //_constraintVector.segment(0,4)         = _controlConstraintVector;
                _constraintVector.segment(0,numRows) = _obstacleConstraintVector;
                
                // Solve the control
                
                Eigen::Vector2d du =del_u;                                                    // We want to solve for this                                        
                try
                {
                    du = QPSolver<double>::solve(d2Ldu2,
                                                 dLdu + d2Ldudx * dx,
                                                 _constraintMatrix,
                                                 _constraintVector,
                                                 du);
                }
                catch (const std::exception &exception)
                {
                    std::cout<<"\n===========ERRORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR================\n";
                    return du;

                    throw std::runtime_error(std::string(exception.what()) + " "
                                            "Failed on recursion no. " + std::to_string(i) + " "
                                           "at step no. " + std::to_string(j) + ".");

                    
                }
                std::cout<<"\n ====== Change in Velocity========== "<<du(0)<<"\t"<<du(1);
                _predictedStates[j].velocity += du ;                                                 // Increment control input

                double norm = du.norm();

                if (norm > largestStepChange) largestStepChange = norm;                             // Save the largest
                
                RobotLibrary::Model::Pose2D tempNewPose = RobotLibrary::Model::DifferentialDrive::predicted_pose(_predictedStates[j].pose,
                                                                   _predictedStates[j].velocity,
                                                                   _controlFrequency);
                
                Eigen::Vector3d new_E  = tempNewPose.error(desiredStates[j+1].pose);
                std::cout<<"\n========New Error ========== \n"<<new_E;
                //dfdx = configuration_jacobian(currentPose, _predictedStates[j].velocity  , _controlFrequency); 
                costateVector = -new_E.transpose() * K *dfdx  + costateVector.transpose()*dfdx; 
            }
        }
        
        // Forward recursions
        for (int j = 0; j < _predictionSteps; ++j)
        {       
            _predictedStates[j+1].pose =
            RobotLibrary::Model::DifferentialDrive::predicted_pose(_predictedStates[j].pose,
                                                                   _predictedStates[j].velocity,
                                                                   _controlFrequency);
                                              
            _predictedStates[j+1].covariance =
            RobotLibrary::Model::DifferentialDrive::predicted_covariance(_predictedStates[j].pose,
                                                                         _predictedStates[j].velocity,
                                                                         _predictedStates[j].covariance,
                                                                         _controlFrequency);
        }
        
        if (largestStepChange < _threshold) break;                                                  // Break early if step change is tiny
    }
    
    return _predictedStates[0].velocity;                                                            // Only return the 1sts
}

} } // namespace
