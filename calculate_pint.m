function pint=calculate_pint(W)

W_dimension = size(W);
n_sections =W_dimension(1); %the number of cross sections along this river segment
b_loga = zeros(n_sections, 2); %initiate a matrix to store rating coefficients
    
mean_log_Q_sections = zeros(n_sections,1);
mean_log_W_sections = zeros(n_sections,1);


% Loop through each cross section:
%figure;
for i_section = 1:n_sections
   
    W_this_section = W(i_section, :)';
    %hold on; plot(log10(Q_this_section), log10(W_this_section))
    
    %re order W_this_section
    [W_this_section, sort_ind] = sort(W_this_section);
  
    
    %Re-scale log(Q) to be [0, 10]
    fQ_min = 0;
    fQ_max = 10;
    log_W_this_section = log10(W_this_section);
    fW_min = mean(log_W_this_section(1:2));
    fW_max = mean(log_W_this_section((length(log_W_this_section)-1):length(log_W_this_section))); 
     fQ = (log_W_this_section-fW_min).*(fQ_max-fQ_min)./(fW_max-fW_min); 
    %figure; plot(log_W_this_section);
    %hold on; plot(fQ, log_W_this_section)
    
    %fit in a and b: logw = loga + blogQ

    %Rescaled pseudo AHG rating curves:
    X = [ones(size(fQ)) fQ];
    fake_coef = polyfit(fQ,log_W_this_section,1);
    fake_b_loga(i_section,:)= fake_coef; % store [b, log10(a)] for this cross-sectional rating


    mean_log_W_sections(i_section) = mean(log10(W_this_section));
    

end


fake_b = fake_b_loga(:,1); fake_loga = fake_b_loga(:,2); 

% Calculate all rating intersections for this river segment

%Rescaled pseudo AHG rating curves
n_line = size(fake_b_loga);
fake_XY_intersections = [];
for i_line = 1:n_line(1)
    for j_line = 1:n_line(1)
        if j_line>i_line
            p1 = [fake_b(i_line), fake_loga(i_line)]; p2 = [fake_b(j_line), fake_loga(j_line)];
            %calculate this intersection
            x_intersect = fzero(@(x) polyval(p1-p2,x),3);
            y_intersect = polyval(p1,x_intersect);
            fake_XY_intersections=[fake_XY_intersections, [x_intersect;y_intersect]];
            %hold on; plot(x_intersect, y_intersect, 'rO', 'MarkerSize',2)
        end
    end
end
fake_logQ_intersected = fake_XY_intersections(1,:); fake_logW_intersected = fake_XY_intersections(2,:);

%
%check the number of intersections within observation ranges:

[fake_within_range_indices]= ...
    find( (fake_logQ_intersected >= fQ_min & fake_logQ_intersected <= fQ_max) & (fake_logW_intersected > log10(min(min(W))) & fake_logW_intersected < log10(max(max(W)))) );
pint = length(fake_within_range_indices)/length(fake_logQ_intersected);





