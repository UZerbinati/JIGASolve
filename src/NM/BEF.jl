function BernstainBasisDerivatives(p,q,s,t)
	#####################################################################
	#
	# Input: p,q degree of 2D Bernstain base function
	# 	 s,t point of IxI
	# Output: value of the Bernstain bases function at (s,t)
	# 	  value of the Bernstain bases function derivatives at (s,t)
	# Warning: we are working with modified parent Bernstain basis,
	# 	   use PAT 
	#
	####################################################################
	B = [];
	∂B = zeros((p+1)*(q+1),2);
	h = 0.001;
	k=1;
	for i in 0:p
		for j in 0:q
			P = zeros(p+1)
			P[i+1] = 1
			CS = BezierCurve(p,P);
			Q = zeros(q);
			Q = zeros(q+1)
			Q[j+1] = 1;
			CT = BezierCurve(q,Q);
			#Computing using DeCasteljau the value of the Bernstain base
			B = push!(B,DeCasteljau(CS,PAT(s))*DeCasteljau(CT,PAT(t)));
			#Computing the derivatives using finite difference
			∂B[k,1] = (1/h)*(DeCasteljau(CS,PAT(s+h))*DeCasteljau(CT,PAT(t))-DeCasteljau(CS,PAT(s))*DeCasteljau(CT,PAT(t)));
			∂B[k,2] = (1/h)*(DeCasteljau(CS,PAT(s))*DeCasteljau(CT,PAT(t+h))-DeCasteljau(CS,PAT(s))*DeCasteljau(CT,PAT(t)));
			k=k+1;
		end 
	end
	return B, ∂B;
end
function BEFShape2D(P,e,Ω,Q,CC,WW,WWb,Nen)
	#####################################################################
	#
	# Input:  P=(ξ,η), quadrature point
	# 	  e, local bezier element number
	# 	  Ω, Bezier Surface
	# 	  Q, polynomial degree
	#	  CC, Bezier extraction operator
	#	  WWb, weights of the Bezier surface
	#	  WW, weights of the NURBS surface
	#	  Nen, number of shape functions
	# Output: shape function values on the quadrature point
	# 	  shape function derivatives value on the quadrature point
	# 	  Jacobian of the shape function trasformation
	# Warning: we are working with modified parent Bernstain basis,
	# 	   use PAT 
	#
	####################################################################
	p = Q[1]; q=Q[2];
	NCPT = (p+1)*(q+1);
	#println(NCPT)
	R = zeros(Nen);
	wb = 0;
	dwb_dξ = zeros(2);
	dR_dξ = zeros(Nen,2);
	dR_dx = zeros(Nen,2);
	dx_dξ = zeros(2,2);
	dξ_dx = zeros(2,2);
	J = 0;
	ξ = PAT(P[1]); #Quadrature point
	η = PAT(P[2]); #Quadrature point
	BΩ = reshape(Ω.B',3,(Ω.p+1)*(Ω.q+1));
	B,∂B = BernstainBasisDerivatives(Ω.p,Ω.q,P[1],P[2]); #Precomputed shape function and derivatives for parent domain.
	#Initializing local operator for the global one that we have as argument
	C = CC[e];
	#Use Bernstein basis to compute the weight functions.
	for a in 1:NCPT
		wb = wb+ B[a]*(diag(WWb[e])[a]);
		dwb_dξ[1] = dwb_dξ[1] + ∂B[a,1]*(diag(WWb[e])[a]);
		dwb_dξ[2] = dwb_dξ[2] + ∂B[a,2]*(diag(WWb[e])[a]);
	end
	IENM = IEN(p,q);
	#Using Re(ξ) = WeCe (Be(ξ)/Wb(ξ)) to compute element shape functions and derivatives with respect to the parent domain.
	for a in 1:Nen
		for b in 1:NCPT
			R[a]=R[a]+WW[Int(IENM[a,e])]*C[a,b]*(B[b]/wb);
			for i in 1:2
				dR_dξ[a,i] = dR_dξ[a,i] + WW[Int(IENM[a,e])]*C[a,b]*(∂B[b,i]/wb-(dwb_dξ[i]*B[b])/(wb*wb));
			end
		end
	end
	#Compute the derivative of the mapping from the parent domain to the physical space.
	for a = 1:NCPT
		for i = 1:2
			for j = 1:2
				dx_dξ[i,j] = dx_dξ[i,j] + BΩ[:,a][i]*∂B[a,j];
			end
		end
	end
	#Computing inverse
	dξ_dx[1,1] = dx_dξ[2,2]; 
	dξ_dx[1,2] = -dx_dξ[1,2];
	dξ_dx[2,1] = -dx_dξ[2,1]; 
	dξ_dx[2,2] = dx_dξ[1,1];
	#println(dξ_dx);
	#println(dx_dξ);
	#println(det(dx_dξ));
	dξ_dx = (1/det(dx_dξ))*dξ_dx;
	
	#Using (∂Re(ξ)/∂xi) = ∑(∂Re(ξ)/∂ξj)(∂ξj/∂xi) to compute the derivatives of the element shape function with respect to the physcial coordinate.
	for a = 1:Nen
		for i = 1:2
			for j = 1:2
				dR_dx[a,i] = dR_dx[a,i] + dR_dξ[a,j]*dξ_dx[j,i];
			end
		end
	end
	J = det(dx_dξ);
	return R,dR_dx,J;
end
