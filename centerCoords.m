function S = centerCoords(S)
% center coordinate system of state such that eyes are at [0, 0]

S.tgt_px_fx = S.tgt_px_fx - S.eye_px_filt_trl;
S.tgt_px_lo = S.tgt_px_lo - S.eye_px_filt_trl;
S.tgt_px_hi = S.tgt_px_hi - S.eye_px_filt_trl;

S.tgt_py_fx = S.tgt_py_fx - S.eye_py_filt_trl;
S.tgt_py_lo = S.tgt_py_lo - S.eye_py_filt_trl;
S.tgt_py_hi = S.tgt_py_hi - S.eye_py_filt_trl;

S.eye_px_filt_trl = S.eye_px_filt_trl - S.eye_px_filt_trl;
S.eye_py_filt_trl = S.eye_py_filt_trl - S.eye_py_filt_trl;

end