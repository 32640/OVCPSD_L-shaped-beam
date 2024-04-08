function [xPhysSns]=HeavisideSns(xTilde,beta,eta)
xPhysSns=-(beta.*(tanh(beta.*(eta - xTilde)).^2 - 1))./(tanh(beta.*eta) - tanh(beta.*(eta - 1)));