function [scalar, vector] = quatmult(s1, v1, s2, v2)
scalar = s1*s2 - dot(v1, v2);
vector = s1*v2 + s2*v1 + cross(v1,v2);