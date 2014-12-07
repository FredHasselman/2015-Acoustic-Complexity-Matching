%%  Keep vector graphics support of 'Renderer', 'painters'
% This will change all patch objects in figure_h to an indexed colormap containing all the unique colors

function patch2ind(figure_h)
objects = findall(figure_h,'Type','patch');
cmapID=-1*ones(1,3);
for o = 1:numel(objects)
 if size(get(objects(o),'EdgeColor'),2)==3
  cmapID = [cmapID;get(objects(o),'EdgeColor')];
 elseif get(objects(o),'MarkerEdgeColor')=='auto'
  cmapID = [cmapID;get(get(get(objects(o),'Parent'),'Parent'),'ColorOrder')];
 else
  cmapID = [cmapID;get(figure_h,'ColorMap')];
 end
end
cmapID = [cmapID; get(figure_h,'ColorMap')];
[Y,newmap] = cmunique(cmapID(cmapID>-1));
ind = 1:size(newmap,1);

for o = 1:numel(objects)
 fv = get(objects(o),'FaceVertexCData');
 if size(fv,2)==3
  if any(ismember(newmap,fv,'rows'))
   set(objects(o),'CDataMapping','direct','FaceVertexCData',ind(ismember(newmap,fv,'rows')))
  end
  
  if size(get(objects(o),'EdgeColor'),2)==3
   ec = get(objects(o),'EdgeColor');
   if any(ismember(newmap,ec,'rows'))
    set(objects(o),'CDataMapping','direct','FaceVertexCData',ind(ismember(newmap,ec,'rows')))
   end
  end
 end
 set(figure_h,'ColorMap',newmap,'Renderer','painters');
end
