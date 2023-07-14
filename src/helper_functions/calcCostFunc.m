function cost_out = calcCostFunc(ave_Consistency_in,n)

%{

    Author: Nicholas Ott

    inputs:
        1) ave_NEES_in -> ave Consistency value for N montecarlo runs over T
        discrete time steps

        2) n -> system order (size of state vector) or size of measurement
        vector

    equation can be found in section B Stochastic Costs for
    Consistency-Based Filter Auto-Tuning in paper

%}


sum_Consistency = sum(ave_Consistency_in);
logarithm = log10(sum_Consistency/n);

cost_out = sqrt(logarithm^2);


end