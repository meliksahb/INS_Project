# INS_Project
INS project using ETH Zurich Micro Aerial Data


# First plot: 3D Trajectories of Each Estimation and Ground Truth
This plot is wrong because it cannot plot other trajectory estimations than pure mechanization eqn solution and ground truth plotting is wrong.

## This plot should be fix!!

![1](https://github.com/user-attachments/assets/ca2a2a23-f4f2-47f2-b4cc-f5da1c43c158)

# Second plot: Horizontal Trajectory Comparison
This plot is also faulty because it cannot plot other trajectory estimations and ground truth than  pure mechanization eqn

## This plot should be fix!!

![2](https://github.com/user-attachments/assets/ed82eac8-3d0c-43dc-8f6e-922c87341f7a)

# Third plot: EKF & UKF Position
EKF and UKF solutions is not in ±3σ boundries for east and north. This should be fixed.

## This plot should be fix!!

![3](https://github.com/user-attachments/assets/e1509452-a261-4668-b587-ed8d24de7cc1)

# Fourth plot: EKF & UKF Velocity
EKF and UKF solutions is in the ±3σ boundries and seems logical.

![4](https://github.com/user-attachments/assets/9478ce48-cf13-4bf6-a012-3551f70cf3ed)

# Fifth plot: EKF & UKF & Mechanization Comparison
I don't know if it is correct or not. But I think UKF solution should be more accurate than EKF. In this plot EKF solution is more accurate. So it should be checked.

![5](https://github.com/user-attachments/assets/fdad028c-238c-4cb5-9ecc-6713afff9aa5)

# Sixth plot: Mechanization Error
This error graph seem correct.

![6](https://github.com/user-attachments/assets/630594dc-aba0-4206-9f78-5679ce1a957a)

# Seventh plot: Loosely vs Tightly Coupled for EKF
Graphs seems incorrect and hroizontal trajectory plot cannot plot gorund truth and Loosely Coupled EKF.

![7](https://github.com/user-attachments/assets/242384d2-702c-411a-8dac-0d85940177ad)

# Eigth plot: Accel Bias Estimation
I estimated biases like this. I didn't check if they are shared model of accel and gyros.

![8](https://github.com/user-attachments/assets/3b2d80d3-5e99-4451-b52d-2073da6ccad6)

# Ninth plot: EKF & UKF NIS Comparison
Normalized Innovation Squared plots.

![9](https://github.com/user-attachments/assets/1b822e0e-da4c-4146-bc14-86874826f6b7)
