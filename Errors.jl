# Compute energy error, phase-space error, and position error between numerical and exact solutions
function errors(H::Function, u0, uNum, uExt)
    diff_u = uNum - uExt
    H_u = mapslices(H, uNum, dims=1)
    err_H = H_u' .- H(u0)

    err_qp = mapslices(norm, diff_u, dims=1)
    err_q = mapslices(norm, diff_u[1:2, :], dims=1)
    return err_qp, err_q, err_H
end
