#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <cmath>

using namespace cv;

/**
 * @brief Core image-processing function for plate bending angle estimation.
 * 
 * @param[in]  inputFrame    Raw BGR frame captured from the camera system.
 * @param[in]  loopCounter   Current iteration index of the vision processing loop.
 * @param[in,out] ioAngleInitial Reference to the baseline initial angle (degrees). 
 *                               Calculated and updated when loopCounter == 15.
 * @param[out] outBendAngle  Calculated instantaneous bending angle (degrees) relative to baseline.
 * @param[in]  debugMode     Optional flag to display OpenCV debug visualization windows.
 * @return true if processing succeeded, false if frame was invalid or points could not be fit.
 */
bool image_processor(const cv::Mat& inputFrame, 
                         int loopCounter, 
                         double& ioAngleInitial, 
                         double& outBendAngle,
                         bool debugMode = false) 
{
    if (inputFrame.empty()) {
        std::cerr << "[Vision Error] Input frame is empty." << std::endl;
        return false;
    }

    try {
        // --- 1. Initial Baseline Angle Calibration (Loop 15) ---
        if (loopCounter == 15) {
            cv::Rect cropRectInitial(50, 200, 250, 320);
            cv::Rect imgBounds(0, 0, inputFrame.cols, inputFrame.rows);
            cv::Rect safeInitial = cropRectInitial & imgBounds;

            if (safeInitial.area() > 0) {
                cv::Mat croppedInitial = inputFrame(safeInitial).clone();
                cv::Mat whiteMaskInitial;
                cv::inRange(croppedInitial, cv::Scalar(200, 200, 200), cv::Scalar(255, 255, 255), whiteMaskInitial);

                std::vector<cv::Point> pointsInitial;
                cv::findNonZero(whiteMaskInitial, pointsInitial);

                if (!pointsInitial.empty()) {
                    cv::Vec4f lineParamsInitial;
                    cv::fitLine(pointsInitial, lineParamsInitial, cv::DIST_L2, 0, 0.01, 0.01);
                    
                    double vx = lineParamsInitial[0];
                    double vy = lineParamsInitial[1];
                    ioAngleInitial = -atan2(vy, vx) * 180.0 / CV_PI;
                } else {
                    std::cerr << "[Vision Warning] No edge points found for initial calibration." << std::endl;
                }

                if (debugMode) {
                    cv::imshow("Cropped Initial", croppedInitial);
                    cv::imshow("White Mask Initial", whiteMaskInitial);
                }
            }
        }

        // --- 2. Real-Time Bending Angle Calculation ---
        cv::Rect cropRect(50, 200, 250, 320); 
        cv::Rect imgBounds(0, 0, inputFrame.cols, inputFrame.rows);
        cv::Rect safeRect = cropRect & imgBounds;

        if (safeRect.width == cropRect.width && safeRect.height == cropRect.height) {
            cv::Mat cropped = inputFrame(safeRect).clone();
            
            // Mask white plate surface/edges
            cv::Mat whiteMask;
            cv::inRange(cropped, cv::Scalar(180, 180, 180), cv::Scalar(255, 255, 255), whiteMask);
            
            // Extract edge coordinates & fit line
            std::vector<cv::Point> points;
            cv::findNonZero(whiteMask, points);

            if (!points.empty()) {
                cv::Vec4f line;
                cv::fitLine(points, line, cv::DIST_L2, 0, 0.01, 0.01);
                
                double rawAngle = -atan2(line[1], line[0]) * 180.0 / CV_PI;
                outBendAngle = rawAngle - ioAngleInitial;
            } else {
                return false; // No points detected in ROI
            }

            if (debugMode) {
                cv::imshow("Cropped View", cropped);
                cv::imshow("White Masked", whiteMask);
                cv::waitKey(1); 
            }
        }
        return true;

    } catch (const cv::exception& e) {
        std::cerr << "[OpenCV Error]: " << e.what() << std::endl;
        return false;
    }
}
