function hit=interactPt(refPts,hid,ptNameList,color)
        hit = impoint(gca,refPts(hid,:));
        fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(hit,fcn);
        setColor(hit,color);
%         setString(hit,ptNameList{hid});
%         addNewPositionCallback(hit,@(h) title([ptNameList{hid},sprintf('(%1.0f,%1.0f)',h(1),h(2))]));
        title([ptNameList{hid}]);
end