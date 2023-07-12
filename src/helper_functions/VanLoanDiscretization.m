function Qout = VanLoanDiscretization(Win,Bin,Ain,dt)

%{

    Author: Nicholas Ott

    Function serves to discretize a continuous process covariance matrix
    using Van Loan's Method


    Win -> Continuous Process Noise Uncertainty
    Bin -> Continuous Input Matrix
    Ain -> Continuous STM
    dt  -> Time Step

%}
    [n,~] = size(Ain);

    A = [-Ain Bin*Win*Bin';zeros(n,n) Ain']*dt;
    B = expm(A);

    lower_right_partition = B(n+1:2*n,n+1:2*n);
    upper_right_partition = B(1:n,n+1:2*n);

    Qout = lower_right_partition'*upper_right_partition;

end