function [fig] = Plot_Param_Estimates(param_estimates,model,observations,true_params,y0,tspan,fill_params,obj_fun)
    if isempty(obj_fun)
        obj_fun = @OLS;
    end
    [~,y_true] = ode45(model,tspan,y0,[],fill_params(true_params));
    true_data = observations(y_true);
    fig = figure;
    hold on
    scatter(param_estimates(1,:),param_estimates(2,:),10,'blue',"filled");
    xl = xlim;
    yl = ylim;
    P1_space = linspace(xl(1),xl(2),100);
    P2_space = linspace(yl(1),yl(2),100);
    [P1,P2] = meshgrid(P1_space,P2_space);
    s = size(P1);
    cols = s(2);
    rows = s(1);
    Z = zeros(s);
    for c=1:cols
        for r = 1:rows
            params = fill_params([P1(r,c),P2(r,c)]);
            if any(params < 0)
                err = 1e10;
            else
                [~,y_new] = ode45(model,tspan,y0,[],params);
                data_new = observations(y_new);
                err = obj_fun(true_data,data_new);
            end
            Z(r,c) = err;
        end
    end
    contour(P1,P2,Z,'LineWidth',1.5);
    scatter(true_params(1),true_params(2),500,"red","filled","pentagram");
    ylim(yl);
    hold off
end
function er = OLS(data1,data2)
    er = sum((data1 - data2).^2,'all');
end