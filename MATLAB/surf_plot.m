function h = surf_plot( X, name, ext )
    h=surf(X);
    
    axis tight;
    caxis([-10,10]);
    zlim([-10,10]);
    %set(gca,'visible','off');
    title(name);
    saveas(gcf,['plots/',name,'.',ext]);
end

