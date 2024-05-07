function S = centerCoords(S)
% center coordinate system of state such that eyes are at [0, 0]

S.tgt_px_fx = subtractWrapper(S.tgt_px_fx, S.eye_px_filt_trl);
S.tgt_px_lo = subtractWrapper(S.tgt_px_lo, S.eye_px_filt_trl);
S.tgt_px_hi = subtractWrapper(S.tgt_px_hi, S.eye_px_filt_trl);

S.tgt_py_fx = subtractWrapper(S.tgt_py_fx, S.eye_py_filt_trl);
S.tgt_py_lo = subtractWrapper(S.tgt_py_lo, S.eye_py_filt_trl);
S.tgt_py_hi = subtractWrapper(S.tgt_py_hi, S.eye_py_filt_trl);

S.eye_px_filt_trl = subtractWrapper(S.eye_px_filt_trl, S.eye_px_filt_trl);
S.eye_py_filt_trl = subtractWrapper(S.eye_py_filt_trl, S.eye_py_filt_trl);

    function diff = subtractWrapper(a,b)
        diff = a-b; 
        if isnan(diff)
            % inf-inf 
            diff=0;
        end
    end

end