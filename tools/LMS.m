function delTheta = LMS(theta,x,y)

delTheta = x'*(y - x*theta);

end