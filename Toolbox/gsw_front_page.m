function gsw_front_page

% Front page to the Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 

if exist('showdemo.m','file') == 2
    showdemo gsw_front_page
else
    [html_file] = which('gsw_front_page.html');
    try
        web ([html_file],' -helpbrowser')
    catch
        try
            web ([html_file],' -browser')
        catch
            disp('Enter the following address into your web browser:')
            disp([html_file])
        end
    end
end

end