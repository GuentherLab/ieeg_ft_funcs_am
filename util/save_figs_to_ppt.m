% 1. Find and sort all open figures
allFigs = findall(groot, 'Type', 'figure');
if isempty(allFigs)
disp('No open figures found to export.');
return;
end
[~, sortIdx] = sort([allFigs.Number]);
allFigs = allFigs(sortIdx);
% 2. Open a connection to PowerPoint via ActiveX
try
pptApp = actxserver('PowerPoint.Application');
pptApp.Visible = 1; % Make PowerPoint visible to track progress
catch ME
error('Could not start PowerPoint. Ensure PowerPoint is installed and you are running Windows.');
end
% 3. Create a new presentation
presentations = pptApp.Presentations;
presentation = invoke(presentations, 'Add');
% Temporary path to cycle image exports
tempImg = fullfile(pwd, 'temp_fig_export.png');
% 4. Loop through figures and inject into slides
for i = 1:numel(allFigs)
figHandle = allFigs(i);
% Force update the figure and export to an image
drawnow;
print(figHandle, tempImg, '-dpng', '-r150'); % 150 DPI matches screen sizing nicely
% Add a blank slide using invoke to avoid toolbox conflicts
% 12 = ppLayoutBlank constant in PowerPoint
slideIndex = presentation.Slides.Count + 1;
slide = invoke(presentation.Slides, 'Add', slideIndex, 12);
% Insert the figure image onto the slide
% Parameters: AddPicture(FileName, LinkToFile, SaveWithDocument, Left, Top, Width, Height)
% Set Width/Height to -1 to maintain the image's original aspect ratio
picture = invoke(slide.Shapes, 'AddPicture', tempImg, 'msoFalse', 'msoTrue', 50, 50, -1, -1);
% Set up slide title string
if isempty(figHandle.Name)
slideTitle = sprintf('Figure %d', figHandle.Number);
else
slideTitle = sprintf('Figure %d: %s', figHandle.Number, figHandle.Name);
end
% Add a title text box to the top of the blank slide
% Parameters: AddTextbox(Orientation, Left, Top, Width, Height)
% 1 = msoTextOrientationHorizontal
textBox = invoke(slide.Shapes, 'AddTextbox', 1, 50, 10, 500, 40);
textBox.TextFrame.TextRange.Text = slideTitle;
textBox.TextFrame.TextRange.Font.Size = 24;
textBox.TextFrame.TextRange.Font.Bold = 'msoTrue';
end
% 5. Save and clean up
savePath = fullfile(pwd, 'All_Figures_ActiveX.pptx');
invoke(presentation, 'SaveAs', savePath);
invoke(presentation, 'Close');
invoke(pptApp, 'Quit');
delete(pptApp);
% Delete the temporary image asset
if exist(tempImg, 'file')
delete(tempImg);
end
fprintf('Successfully built and saved: %s\n', savePath);