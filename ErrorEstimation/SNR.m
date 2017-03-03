function [ snr ] = SNR( patmat_tx,patmat_rx )

    snr = var(patmat_tx)/var(patmat_tx - patmat_rx);

end

