figure; imagesc(Qtbls{3}); colorbar;
ylabel('State')
yticks(1:4)
yticklabels({'Choice', 'At Low', 'At High', 'Waiting'})
xlabel('Action')
xticks(1:3)
xticklabels({'Go To Low', 'Go To High', 'Fixate'})
title('Q Table')