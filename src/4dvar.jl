# assume we already have constant defined with
# indices of the measurements
function ObsOpMeas(u, idx_meas)   # <--- observation operator for direct chemical measurements
    return u[idx_meas, :]
end
