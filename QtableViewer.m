%% QtableViewer 
% take a specified Q table from the MDP "decision space" formulation and
% display it with labels.
N = length(Qtbls); 
H = floor(sqrt(N)); W = ceil(N/H); 
figure; 
for n = 1:N
    subplot(H,W,n)
    imagesc(Qtbls{n}); colorbar;
    ylabel('State')
    yticks(1:4)
    yticklabels(StateNames)
    xlabel('Action')
    xticks(1:3)
    xticklabels(ActNames)
    title('Q Table')
end