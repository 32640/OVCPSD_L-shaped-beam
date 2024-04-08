function [xPhys]=Heaviside(xTilde,beta,eta)
xPhys=(tanh(beta.*eta)+tanh(beta.*(xTilde-eta)))./(tanh(beta.*eta)+tanh(beta.*(1-eta)));

