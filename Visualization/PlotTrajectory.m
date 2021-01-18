function PlotTrajectory(trajectory)
% Function that plots a trajectory
%
% Inputs:
%   trajectory = [readout_samples dimensions interleaves] OR
%   trajectory = [dimensions readout_samples interleaves]
%
% Warning:
%   Don't plot too many interleaves (3D, >200), this may take a long time.
%
%
% Niek Huttinga, UMC Utrecht, 2020

if size(trajectory,1)<=3
    warning('Detected BART trajectory input format... permuting');
    trajectory = permute(trajectory,[2 1 3]);
end

if size(trajectory,2)<3
    dim = 2;
    % add z=0 to the 2D trajectory for visualization
    trajectory(:,3,:) = zeros(size(trajectory(:,1,:)));
else
    dim=3;
end


for j=1:4
    subplot(2,2,j)
    for i=1:size(trajectory,3)
        scatter3(trajectory(:,1,i),trajectory(:,2,i),trajectory(:,3,i),'LineWidth',1,'Marker','.');hold on;
    end
    xlabel('Kx'),ylabel('Ky'),zlabel('Kz');
    if j==1
        view(0,90);
    elseif j==2
        view(0,0);
    elseif j==3
        view(90,0);
    elseif j==4
        view(3);
    end
end

mtit([num2str(dim),'D Trajectory']);
grid on;
end
