import numpy as np

def richardson_extrapolation(u,Q,relFrac):
    '''
    Conducts repeated richardson extrapolation on
    the matrix u, Q is the fracion between the step lengths, relFrac is the
    condition on the stability of the fraction between the delta-values
    for starting richardson extrapolation
    '''
    nn=len(u);
    u_hat=np.zeros((nn,nn));
    du_hat=np.zeros((nn,nn));
    du_hat_frac=np.zeros((nn,nn));
    du_hat_frac_frac=np.zeros((nn,nn));
    u_hat[:,0]=u[:,0]
    h_p=np.zeros(nn)
    # The counter acc_row(i) keeps track of what row dT_hat_frac(:,i) is first
    # to fulfill the relFrac requirement
    acc_row=np.zeros(nn); acc_row[0]=1; ii=0
    while ii<nn and acc_row[ii]<nn-2: # accrow(ii)+2:end-1 demands that acc_row is never larger than n-3
        du_hat[:,ii]=np.concatenate((np.zeros(acc_row[ii]+1),np.diff(u_hat[acc_row[ii]:,ii])),0)
        du_hat_frac[:,ii]=np.concatenate((np.zeros(acc_row[ii]+1), du_hat[acc_row[ii]+1:-1,ii]/du_hat[acc_row[ii]+2:,ii], np.zeros(1)),0)
        du_hat_frac_frac[:,ii]=np.concatenate([np.zeros(acc_row[ii]+1), du_hat_frac[acc_row[ii]+1:-2,ii]/du_hat_frac[acc_row[ii]+2:-1,ii], np.zeros(2)],0)
        jj=acc_row[ii]+1; bolCont=True
        while jj<nn and bolCont:
            if abs(1-du_hat_frac_frac[jj,ii]) < relFrac:
                acc_row[ii+1]=jj
                p=round(np.log(du_hat_frac[jj,ii])/np.log(Q))
                h_p[ii]=p
                u_hat[jj:,ii+1]=u_hat[jj:,ii]+du_hat[jj:,ii]/(Q**p-1)
                bolCont=False
            jj += 1
        if bolCont: acc_row[ii+1]=nn
        ii += 1
    u_accept=u_hat[-1,ii-1]
    e_trunc=abs(u_hat[-1,ii-1]-u_hat[-2,ii-1])
    return u_accept, e_trunc, h_p



def test_richard():
    start = 5
    stop  = 11
    nn = stop - start + 1
    I     = np.zeros((nn,1))
    e     = range(start,stop+1)
    for ii in range(0,nn):
        n=2**e[ii]
        x=np.linspace(0,1,n)
        y=np.exp(x)
        I[ii]=(np.dot((y[1:]+y[0:-1]),np.diff(x))/2)
    print I-np.exp(1)+1
    print richardson_extrapolation(I,2,1e-2)

def test_fd():
    from get_derivative_from_finite_difference import *
    start = 2
    stop  = 6
    nn    = stop - start + 1
    dydx  = np.zeros((nn,1))
    e     = range(start,stop+1)
    base  = 2
    for ii in range(0,nn):
        n=base**e[ii]
        x=np.linspace(0,1,n)
        y=np.exp(x)
        dydx[ii]=get_derivative_from_finite_difference(1,x,y,0.5)
    print dydx-np.exp(0.5)
    print richardson_extrapolation(dydx,base,1e-1)[0]-np.exp(0.5)


if __name__=='__main__':
    test_richard()