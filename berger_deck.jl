using Plots
# gr()
using QuadGK

# constants
mπ=0.14; mρ=0.77; mf0=0.99;
mssq = [[mπ^2, mf0^2], [mπ^2, mρ^2]];
gs = [3.0,6.0];
Msq0 = 1.26^2;
sl = [1.3, 0.7];


# coupling-strength of deck
cdf0 = +10; cdρ = 1;

# definition of functions
function ChewMandestam(s,m1sq,m2sq)
    sth = (sqrt(m1sq)+sqrt(m2sq))^2;
    spth = (sqrt(m1sq)-sqrt(m2sq))^2;
    dmth = sqrt(Complex(sth-s));
    dmpth = sqrt(Complex(spth-s));
    lamb = dmth*dmpth;
    val = -1/s*lamb*(log((dmth+dmpth)/2)-log(m1sq*m2sq)/4)+
        ((m1sq == m2sq) ? -1 : (
            (m1sq-m2sq)/s*log(m1sq/m2sq)-
            (m1sq+m2sq)/(m1sq-m2sq)*log(m1sq/m2sq)-1/2));
    return -2*val/π;
end

Msq = Msq0 + real(
        gs[1]^2*ChewMandestam(Msq0,mssq[1][1],mssq[1][2])+
        gs[2]^2*ChewMandestam(Msq0,mssq[2][1],mssq[2][2]))/(16*π)

# some more functions
λ(x,y,z) = x^2+y^2+z^2 - 2*x*y - 2*y*z - 2*z*x
ρ = [x->(x>(mf0+mπ)^2) ? sqrt(λ(x,mf0^2,mπ^2))/x : 0,
     x->(x>(mρ +mπ)^2) ? sqrt(λ(x,mρ^2, mπ^2))/x : 0]

# Omnes function
function Dm(s,Msq,g,mssq)
    CMf = [ChewMandestam(s,mssq[i][1],mssq[i][2]) for i=1:2]/(16*π)
    Δ = Msq-s-g[1]^2*CMf[1]-g[2]^2*CMf[2]
    α = g[1]^2+g[2]^2
    1./Δ*[g[1] -g[2]*(Msq-s-α*CMf[2]); g[2]  g[1]*(Msq-s-α*CMf[1])]
end

# # plot just a1 1260
# plot(s->imag(Dm(s+0.0001im,Msq,gs,mssq))[1,1], (mπ+mρ)^2, 2)
# plot!(s->real(Dm(s+0.0001im,Msq,gs,mssq))[1,1], (mπ+mρ)^2, 2)
# plot(e->abs(Dm(e^2+0.0001im,Msq,gs,mssq)[1,1])^2*sqrt(λ(e^2,mssq[2][1],mssq[2][2])),
#         (mπ+mρ), 2)

# some tests
# ChewMandestam(1.1+0.1im,0.14^2,0.14^2)
# ChewMandestam(1.1-0.1im,0.14^2,0.14^2)
# ChewMandestam(1.1,0.14^2,0.14^2)
# plot(s->real(ChewMandestam(s+0.00001im,0.14^2,0.14^2)), -1, 1.1)
# plot!(s->imag(ChewMandestam(s+0.00001im,0.14^2,0.14^2)), -1, 1.1)
# plot!(s->(s>0.28^2) ? sqrt(1-4*0.14^2/s) : 0, -1, 1.1)
# heatmap(-1:0.03:1,-1:0.03:1,
#         (x,y)->imag(ChewMandestam(x+1im*y,0.14^2,0.14^2)))

# basic functions

# plot inverse Omnes matrix
# plot(s->imag(-inv(Dm(s+0.0001im,Msq,gs,mssq)))[1,1], (mπ+mρ)^2, 2)
# plot!(s->imag(-inv(Dm(s+0.0001im,Msq,gs,mssq))[1,2]), (mπ+mρ)^2, 2)
# plot!(s->imag(-inv(Dm(s+0.0001im,Msq,gs,mssq))[2,1]), (mπ+mρ)^2, 2)
# plot!(s->imag(-inv(Dm(s+0.0001im,Msq,gs,mssq))[2,2]), (mπ+mρ)^2, 2)
bDeck = [s->(s>(mf0+mπ)^2) ? cdf0*10*sqrt(λ(s,mf0^2,mπ^2)/s)*(s-1.39^2)*exp(-sl[1]*s)/s^2 : 0,
         s->(s>(mρ+mπ)^2) ? cdρ*8*exp(-sl[2]*s) : 0]

pl = plot(layout = (2,2), size=(800,800))
plot(pl[1], xlab="\$M_{3\\pi}\\textrm{ (GeV)}\$", ylab = "Intensity (a.u.)")
plot!(pl[1], e->10*bDeck[1](e^2)^2*ρ[1](e^2), (mπ+mf0), 1.8, lab = "f0pi x 10")
plot!(pl[1], e->bDeck[2](e^2)^2*ρ[2](e^2), (mπ+mρ), 1.8, lab = "rho pi")

# plot(e->10*bDeck[1](e^2)^2*ρ[1](e^2), (mπ+mf0), 1.8, lab = "f0pi x 10")

# # f0 pi Deck
# plot(s->imag(bDeck[1](s)), 1, 4)
# plot!(s->real(bDeck[1](s)), 1, 4)

# # rho pi Deck
# plot(s->imag(bDeck[2](s)), 0, 4)
# plot!(s->real(bDeck[2](s)), 0, 4)

# unitarized deck, U^D function
UDeck(s) = 1/π*quadgk(sp->(gs[1]*ρ[1](sp)/(16*π)*bDeck[1](sp)+
                           gs[2]*ρ[2](sp)/(16*π)*bDeck[2](sp))/(sp-s-0.00001im),
                    (mρ+mπ)^2, 10)[1]

# plot!(pl[2], xlab="\$M_{3\\pi}^2\\textrm{ (GeV}^2)\$", ylab = "Amplitude (a.u.)", title="Unitarization function")
# plot!(pl[2], s->imag(UDeck(s)), 0.8, 1.8, lab="imag")
# plot!(pl[2], s->real(UDeck(s)), 0.8, 1.8, lab="real")

# Plot Unitarized deck for two channels
function Fu(s)
        Dmv = Dm(s+0.0001im,Msq,gs,mssq)
        bD = [bDeck[1](s), bDeck[2](s)]
        bD + UDeck(s)*Dmv[:,1]
end

# plot(s->real(Fu(s)[1]),0.8, 1.8)
# plot!(s->imag(Fu(s)[1]),0.8, 1.8)

plot!(pl[3],title="\$1^{++} f_0\\pi\\, P\\textrm{-wave intenities}\$")
plot!(pl[3],e->abs(bDeck[1](e^2))^2*sqrt(λ(e^2,mssq[1][1],mssq[1][2])), (mπ+mf0), 2, lab="Deck")
plot!(pl[3],e->abs(Fu(e^2)[1])^2*sqrt(λ(e^2,mssq[1][1],mssq[1][2])), (mπ+mf0), 2, lab="UniDeck")

# plot(s->real(Fu(s)[2]),0.8, 1.8)
# plot!(s->imag(Fu(s)[2]),0.8, 1.8)
plot!(pl[4],title="\$1^{++} \\rho\\pi\\, S\\textrm{-wave intenities}\$")
plot!(pl[4],e->abs(bDeck[2](e^2))^2*sqrt(λ(e^2,mπ^2,mρ^2)), (mπ+mρ), 2, lab="Deck")
plot!(pl[4],e->abs(Fu(e^2)[2])^2*sqrt(λ(e^2,mπ^2,mρ^2)), (mπ+mρ), 2, lab="UniDeck")

function arg(x)
    y = x*exp(-1im*π)
    atan2(imag(y), real(y))+π
end

# phase difference
plot!(pl[2],e->180/π*arg(Fu(e^2)[1]*conj(Fu(e^2)[2])), (mπ+mρ), 2, lab="Phase difference, f0-rho")
plot(pl)
# savefig("/tmp/pm10s1.pdf")
