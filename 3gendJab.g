FF:=Rationals;
alp:=Indeterminate(FF,"alp");;
bet:=Indeterminate(FF,"bet");;
x:=Indeterminate(FF,"x");;
y:=Indeterminate(FF,"y");;
z:=Indeterminate(FF,"z");;
p:=Indeterminate(FF,"p");;

F:=FunctionField(FF, 6);;
o:=One(F);;
a:=[1,0,0,0,0,0,0,0,0]*o;;
b:=[0,1,0,0,0,0,0,0,0]*o;;
c:=[0,0,1,0,0,0,0,0,0]*o;;
ab:=[0,0,0,1,0,0,0,0,0]*o;;
ac:=[0,0,0,0,1,0,0,0,0]*o;;
bc:=[0,0,0,0,0,1,0,0,0]*o;;
abc:=[0,0,0,0,0,0,1,0,0]*o;;
bac:=[0,0,0,0,0,0,0,1,0]*o;;
cab:=[0,0,0,0,0,0,0,0,1]*o;;

##Basic epsilon values
epsab:=-(2*x*(alp-1)+(alp+bet))/(alp-bet);;
epsac:=-(2*y*(alp-1)+(alp+bet))/(alp-bet);;
epsbc:=-(2*z*(alp-1)+(alp+bet))/(alp-bet);;

#Filling the gram matrix and adding additional eps, when needed
gram:=[];
Add(gram,o*[1,x,y,x,y,p,p,epsab*y+z-epsab*(p-epsab*(y-z)),epsac*x+z-epsac*(p-epsab*(y-z))]);;
Add(gram,o*[0,1,z,x-epsab*(x-1),p-epsab*(y-z),z,y-epsab*(p-z),p*epsab*(y-z),y-epsab*(p-z)+epsbc*(x-epsab*(x-1)-(p*epsab*(y-z)))]);;
Add(gram,o*[0,0,1,p-epsab*(y-z),y-epsac*(y-1),z-epsbc*(z-1),epsbc*y+x-epsbc*(p-epsbc*(y-x))-epsac*(p-(z-epsbc*(z-1))),epsbc*y+x-epsbc*(p-epsbc*(y-x))-epsac*(p-(z-epsbc*(z-1))),p-epsab*(y-z)]);;
epstacb:=-(2*gram[2][5]*(alp-1)+(alp+bet))/(alp-bet);; #[a]c,b
Add(gram,o*[0,0,0,1,z,y-epsab*(p-z),z,epsab*(gram[2][5]-gram[1][8])+y,epstacb*(1-gram[2][4])+gram[2][8]]);;
epstabc:=-(2*gram[3][4]*(alp-1)+(alp+bet))/(alp-bet);; #[a]b,c
Add(gram,o*[0,0,0,0,1,gram[3][7],gram[3][6],epstabc*(1-gram[3][4])+gram[3][4],epstacb*(z-gram[3][5])+gram[3][8]]);;
epstbac:=-(2*p*(alp-1)+(alp+bet))/(alp-bet);; #[b]a,c=a,[b]c
Add(gram,o*[0,0,0,0,0,1,epstbac*(1-p)+p,gram[3][5],epsbc*(gram[3][4]-gram[2][9])+gram[2][4]]);;
Add(gram,o*[0,0,0,0,0,0,1,epstabc*(gram[3][6]-gram[4][6])+gram[6][9],epstacb*(gram[2][6]-gram[5][6])+gram[6][8]]);;
epstcba:=-(2*(p+epsbc*(x-y))*(alp-1)+(alp+bet))/(alp-bet);; #[c]b,a=b[c]+epsbc(b-c),a
Add(gram,o*[0,0,0,0,0,0,0,1,epstcba*(x-gram[2][9])+z+epsbc*(1-z)+epsac*(gram[6][9]-(gram[4][9]+epsab*(gram[1][9]-gram[2][9])))]);;
Add(gram,o*[0,0,0,0,0,0,0,0,1]);;
epstbca:=epstbac;;
epstcab:=epstcba;;

mas:=[];

#Words of length 2
ba:=ab+epsab*(a-b);;
cb:=bc+epsbc*(b-c);;
ca:=ac+epsac*(a-c);;

#Words of length 3
acb:=abc+epsbc*(ab-ac);;
cba:=cab+epsab*(ca-cb);;
bca:=bac+epsac*(ba-bc);;
ba_b:=epsab*(b-ba)+a;;
ab_a:=epsab*(a-ab)+b;;
bc_b:=epsbc*(b-bc)+c;;
ca_c:=epsac*(c-ca)+a;;
cb_c:=epsbc*(c-cb)+b;;
ac_a:=epsac*(a-ac)+c;;

#Words of length 4

aba_c:=epstabc*(c-ab)+cab;;
aca_b:=epstacb*(b-ac)+bac;;
bab_c:=epstbac*(c-ba)+cba;;
bcb_a:=epstbca*(a-bc)+abc;;
cac_b:=epstcab*(b-ca)+bca;;
cbc_a:=epstcba*(a-cb)+acb;;
acb_a:=aca_b+epsab*(ac_a-acb);;
abc_a:=aba_c+epsac*(ab_a-abc);;
bca_b:=bcb_a+epsab*(bc_b-bca);;
bac_b:=bab_c+epsbc*(ba_b-bac);;
cab_c:=cac_b+epsbc*(ca_c-cab);;
cba_c:=cbc_a+epsac*(cb_c-cba);;

b_ca_c:=epsac*(bc-bca)+ba;;
b_ab_a:=epsab*(ba-ba_b)+b;;
c_ba_b:=epsab*(cb-cba)+ca;;
a_cb_c:=epsbc*(ac-acb)+ab;;
a_bc_b:=epsbc*(ab-abc)+ac;;
a_ca_c:=epsac*(ac-ac_a)+a;;
c_bc_b:=epsbc*(cb-cb_c)+c;;
b_ac_a:=epsac*(ba-bac)+bc;;
c_ab_a:=epsab*(ca-cab)+cb;;
b_cb_c:=epsbc*(bc-bc_b)+b;;
a_ba_b:=epsab*(ab-ab_a)+a;;

Add(mas,[a,((alp-bet)*ab-2*(alp-1)*x*a+(alp+bet)*b)/2,((alp-bet)*ac-2*(alp-1)*y*a+(alp+bet)*c)/2,((alp-bet)*b-2*(alp-1)*x*a+(alp+bet)*ab)/2,((alp-bet)*c-2*(alp-1)*y*a+(alp+bet)*ac)/2,((alp-bet)*abc-2*(alp-1)*p*a+(alp+bet)*bc)/2,((alp-bet)*bc-2*(alp-1)*p*a+(alp+bet)*abc)/2,((alp-bet)*aba_c-2*(alp-1)*gram[1][8]*a+(alp+bet)*bac)/2,((alp-bet)*aca_b-2*(alp-1)*gram[1][9]*a+(alp+bet)*cab)/2]);;

Add(mas,[mas[1][2], b, ((alp-bet)*bc-2*(alp-1)*z*b+(alp+bet)*c)/2, ((alp-bet)*ba_b-2*(alp-1)*gram[2][4]*b+(alp+bet)*ab)/2, ((alp-bet)*bac-2*(alp-1)*gram[2][5]*b+(alp+bet)*ac)/2, ((alp-bet)*c-2*(alp-1)*gram[2][6]*b+(alp+bet)*bc)/2, ((alp-bet)*bab_c-2*(alp-1)*gram[2][7]*b+(alp+bet)*abc)/2, ((alp-bet)*ac-2*(alp-1)*gram[2][8]*b+(alp+bet)*bac)/2, ((alp-bet)*bca_b-2*(alp-1)*gram[2][9]*b+(alp+bet)*cab)/2]);;

Add(mas,[mas[1][3],mas[2][3],c,((alp-bet)*cab-2*(alp-1)*gram[3][4]*c+(alp+bet)*ab)/2,((alp-bet)*ca_c-2*(alp-1)*gram[3][5]*c+(alp+bet)*ac)/2,((alp-bet)*cb_c-2*(alp-1)*gram[3][6]*c+(alp+bet)*bc)/2,((alp-bet)*cab_c-2*(alp-1)*gram[3][7]*c+(alp+bet)*abc)/2,((alp-bet)*cba_c-2*(alp-1)*gram[3][8]*c+(alp+bet)*bac)/2, ((alp-bet)*ab-2*(alp-1)*gram[3][9]*c+(alp+bet)*cab)/2]);;

#Words of length 5
a_bab_c:=epstbac*(ac-ab_a)+acb_a;;
b_aba_c:=epstabc*(bc-ba_b)+bca_b;;
b_aca_b:=epstacb*(b-bac)+ac;;
b_cbc_a:=epstcba*(ba-bc_b)+bac_b;;
b_cba_c:=b_cbc_a+epsac*(b_cb_c-bcb_a);;
a_bab_c:=epstbac*(ac-ab_a)+acb_a;;
a_cac_b:=epstcab*(ab-ac_a)+abc_a;;
a_cbc_a:=epstcba*(a-acb)+cb;;
a_bcb_a:=epstbca*(a-abc)+bc;;
a_bca_b:=bcb_a+epsab*(a_bc_b-abc_a);;
c_aba_c:=epstabc*(c-cab)+ab;;
c_aca_b:=epstacb*(cb-ca_c)+cba_c;;
c_bab_c:=epstbac*(c-cba)+ba;;
c_bcb_a:=epstbca*(ca-cb_c)+cab_c;;
a_cab_c:=a_cac_b+epsbc*(a_ca_c-aca_b);;
a_cba_c:=a_cbc_a+epsac*(a_cb_c-acb_a);;
a_bac_b:=a_bab_c+epsbc*(a_ba_b-aba_c);;
b_cab_c:=(epstcab*(b-bca)+ca)+epsbc*((epsac*(bc-bca)+ba)-bca_b);;
c_bca_b:=(epstbca*(ca-cb_c)+cab_c)+epsab*(c_bc_b-(cba_c+epsac*(cba-cb_c)));;
b_abc_a:=b_aba_c+epsac*(b_ab_a-bab_c);;
c_abc_a:=c_aba_c+epsac*(c_ab_a-cab_c);;

#Words of length 6
a_b_aba_c:=epstabc*(abc-(epsab*(ab-ab_a)+a))+(epstbca*(a-abc)+bc)+epsab*(epsbc*(ab-abc)+ac-abc_a);;
a_b_aca_b:=epstacb*(ab-aba_c)+c;;
a_c_aba_c:=epstabc*(ac-aca_b)+b;;
a_b_cbc_a:=epstcba*(ab_a-a_bc_b)+a_bac_b;;
c_b_aba_c:=epstabc*(cb_c-c_ba_b)+c_bca_b;;
c_b_aca_b:=epstacb*(cb-cba_c)+ca_c;;
b_c_bc_b:=epsbc*(bc_b-b_cb_c)+bc;;
a_c_aca_b:=epstacb*(acb-(epsac*(ac-ac_a)+a))+(epstcba*(a-acb)+cb)+epsac*((epsbc*(ac-acb)+ab)-(aca_b+epsab*(ac_a-acb)));;
b_c_bab_c:=epstbac*(bc-(bca_b+epsab*(bca-bc_b)))+a;;
b_c_bca_b:=epstbca*(bca-b_cb_c)+b_cab_c+epsab*(b_c_bc_b-(b_cba_c+epsac*(bcb_a-b_cb_c)));;
###############

#Words required to derive [a,b]c*[b,a]c.


b_c_b_aba_c:=epstabc*(b_cb_c-(epsab*(bc_b-bcb_a)+bca))+b_c_bca_b;;

a_b_cb_c:=epsbc*(abc-a_bc_b)+ab;;
a_b_cab_c:=(epstcab*(ab-abc_a)+ac_a)+epsbc*((epsac*(abc-abc_a)+ab_a)-a_bca_b);;

a_b_c_bca_b:=epstbca*(abc_a-a_b_cb_c)+a_b_cab_c+epsab*((epsbc*(a_bc_b-a_b_cb_c)+abc)-(a_cba_c+epsac*(acb_a-a_b_cb_c)));;

a_b_c_b_aba_c:=epstabc*(a_b_cb_c-(epsab*(a_bc_b-a_bcb_a)+abc_a))+a_b_c_bca_b;;


#Words required to derive [a,b]c*[c,a]b
a_b_ca_c:=epsac*(abc-abc_a)+ab_a;;
a_b_cba_c:=a_b_cbc_a+epsac*(a_b_cb_c-a_bcb_a);;

b_c_b_aca_b:=epstacb*(bc_b-b_cba_c)+b_ca_c;;
a_b_c_b_aca_b:=epstacb*(a_bc_b-a_b_cba_c)+a_b_ca_c;;

#Words required to derive [b,a]c*[c,a]b

c_a_bc_b:=epsbc*(cab-cab_c)+ca_c;;

c_a_bca_b:=c_bcb_a+epsab*(c_a_bc_b-c_abc_a);;

a_c_bcb_a:=epstbca*(ac_a-a_cb_c)+a_cab_c;;

a_c_a_bc_b:=epsbc*(aca_b-a_cab_c)+a_ca_c;;

a_c_ab_a:=epsab*(ac_a-aca_b)+acb;;

a_c_abc_a:=a_c_aba_c+epsac*(a_c_ab_a-a_cab_c);;

a_c_a_bca_b:=a_c_bcb_a+epsab*(a_c_a_bc_b-a_c_abc_a);;

b_a_cb_c:=epsbc*(bac-bac_b)+ba_b;;


##The final step

b_a_cac_b:=epstcab*(ba_b-b_ac_a)+b_abc_a;;

b_a_ca_c:=epsac*(bac-b_ac_a)+ba;;
b_a_c_aba_c:=epstabc*(bac-b_aca_b)+b;;
b_a_c_ab_a:=epsab*(b_ac_a-b_aca_b)+bac_b;;


b_a_cab_c:=b_a_cac_b+epsbc*(b_a_ca_c-b_aca_b);;

b_a_c_bcb_a:=epstbca*(b_ac_a-b_a_cb_c)+b_a_cab_c;;
b_a_c_a_bc_b:=epsbc*(b_aca_b-b_a_cab_c)+b_a_ca_c;;
b_a_c_abc_a:=b_a_c_aba_c+epsac*(b_a_c_ab_a-b_a_cab_c);;

b_a_c_a_bca_b:=b_a_c_bcb_a+epsab*(b_a_c_a_bc_b-b_a_c_abc_a);;

#Finisihing the multiplication matrix
Add(mas,[mas[1][4],mas[2][4],mas[3][4],ab,((alp-bet)*abc-2*(alp-1)*gram[4][5]*ab+(alp+bet)*ac)/2,((alp-bet)*a_bab_c-2*(alp-1)*gram[4][6]*ab+(alp+bet)*bc)/2,((alp-bet)*ac-2*(alp-1)*gram[4][7]*ab+(alp+bet)*abc)/2, ((alp-bet)*a_b_aba_c-2*(alp-1)*gram[4][8]*ab+(alp+bet)*bac)/2,((alp-bet)*a_b_aca_b-2*(alp-1)*gram[4][9]*ab+(alp+bet)*cab)/2]);;

Add(mas,[mas[1][5],mas[2][5],mas[3][5],mas[4][5],ac,((alp-bet)*a_cab_c-2*(alp-1)*gram[5][6]*ac+(alp+bet)*bc)/2, ((alp-bet)*a_cb_c-2*(alp-1)*gram[5][7]*ac+(alp+bet)*abc)/2, ((alp-bet)*a_c_aba_c-2*(alp-1)*gram[5][8]*ac+(alp+bet)*bac)/2,((alp-bet)*a_c_aca_b-2*(alp-1)*gram[5][9]*ac+(alp+bet)*cab)/2]);;

Add(mas,[mas[1][6],mas[2][6],mas[3][6],mas[4][6],mas[5][6],bc,((alp-bet)*b_c_bab_c-2*(alp-1)*gram[6][7]*bc+(alp+bet)*abc)/2, ((alp-bet)*b_ca_c-2*(alp-1)*gram[6][8]*bc+(alp+bet)*bac)/2, ((alp-bet)*b_c_bca_b-2*(alp-1)*gram[6][9]*bc+(alp+bet)*cab)/2]);;

Add(mas,[mas[1][7],mas[2][7],mas[3][7],mas[4][7],mas[5][7],mas[6][7],abc, ((alp-bet)*a_b_c_b_aba_c-2*(alp-1)*gram[7][8]*abc+(alp+bet)*bac)/2, ((alp-bet)*a_b_c_b_aca_b-2*(alp-1)*gram[7][9]*abc+(alp+bet)*cab)/2]);;

Add(mas,[mas[1][8],mas[2][8],mas[3][8],mas[4][8],mas[5][8],mas[6][8],mas[7][8],bac,((alp-bet)*b_a_c_a_bca_b-2*(alp-1)*gram[8][9]*bac+(alp+bet)*cab)/2]);;

Add(mas,[mas[1][9],mas[2][9],mas[3][9],mas[4][9],mas[5][9],mas[6][9],mas[7][9],mas[8][9],bac]);;

#The product function
prod := function(u,v)
  local i,j,k, ans;
  ans:=o* 0*a;
  for i in [1..9] do
  for j in [1..9] do
  for k in [1..9] do
    ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
  od; od; od;
  return ans;
end;

##### verifying the spectra of ad_a, ad_b, ad_c
ad_a:=o*[mas[1][1], mas[1][2], mas[1][3], mas[1][4], mas[1][5], mas[1][6], mas[1][7], mas[1][8],mas[1][9]];;
 (ad_a-ad_a^0)*(ad_a-alp*ad_a^0)*(ad_a-bet*ad_a^0);
ad_b:=o*[mas[2][1], mas[2][2], mas[2][3], mas[2][4], mas[2][5], mas[2][6], mas[2][7], mas[2][8],mas[2][9]];;
(ad_b-ad_b^0)*(ad_b-alp*ad_b^0)*(ad_b-bet*ad_b^0);
ad_c:=o*[mas[3][1], mas[3][2], mas[3][3], mas[3][4], mas[3][5], mas[3][6], mas[3][7], mas[3][8],mas[3][9]];;
(ad_c-ad_c^0)*(ad_c-alp*ad_c^0)*(ad_c-bet*ad_c^0);


