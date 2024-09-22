function errors(H::Function,u0,uNum,uExt)
    diff_u = uNum-uExt; # matrix contains vectors whose elements are u = (q1,q2,p1,p2) 
       H_u = mapslices(H,uNum,dims=1); # eval H at each vector whose elements are u = (q1,q2,p1,p2) 
     err_H = H_u' .- H(u0);

    err_qp = mapslices(norm,diff_u,dims=1);  # norm eachcolumn of diff_u
    err_q  = mapslices(norm,diff_u[1:2,:],dims=1);
    return err_qp, err_q, err_H; 
end