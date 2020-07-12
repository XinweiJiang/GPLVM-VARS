function [ output_vector, r1 ] = joinvector( vector1, vector2 )

% this function is used to join rows of any two vectors

[r1,c1] = size(vector1);
[r2,c2] = size(vector2);

if c1 == c2
    for i_outer = 1: r1
       output_vector(i_outer,:)= vector1(i_outer,:);
    end
    count = 0;
    for i = r1+1 : r1+r2
        count = count+1;
        output_vector(i,:)= vector2(count,:);
    end
else
    output_vector = [];
end

if r1==0 && r2==0
    output_vector = [];
end

if r1==0 && r2~=0
    output_vector = vector2;
end

if r2==0 && r1~=0
    output_vector = vector1;
end
        

end

