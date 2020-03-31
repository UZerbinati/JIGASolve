S = "((7+3)*(5-2))";
mutable struct BST
	I::Int64;
	O::String;
	L::Union{Missing,BST};
	R::Union{Missing,BST};
end
function TreeEx()
	##################################
	#            0
	#           / \
	#          1   2
	#         / \
	#        3   4
	#             \
	#              5
	###################################       
	T = BST(0,"0",Missing(),Missing())
	N5 = BST(5,"5",Missing(),Missing());
	N3 = BST(3,"3",Missing(),Missing());
	N4 = BST(4,"4",Missing(),N5);
	N1 = BST(1,"1",N3,N4);
	N2 = BST(2,"2",Missing(),Missing());
	T.R = N2; T.L=N1;
	TreePrint([T],5);
end
function TreePrint(S,i)
	Q = []
	if i == 1
		return;
	end
	M = "";
	P=" "^(2^i-1)
	PP=" "^Int((0.5*2^i-1)-1)
	PPP=" "^((2^i-1)-(Int((0.5*2^i-1))))
	Z="-"^Int((0.5*2^i-1)-1)*" "*PPP
	I = [!isa(s,Missing) for s in S]
	if S[I] == []
		return;
	end
	for T in S
		#println(T)
		if isa(T,Missing)
			push!(Q,Missing())
			push!(Q,Missing())
			print("#"*Z)
			M=M*"|"*PP*"\\"*PPP
		elseif !(isa(T.R,Missing)) && !(isa(T.L,Missing))
			push!(Q,T.L)
			push!(Q,T.R)
			print(T.O*Z)
			M=M*"|"*PP*"\\"*PPP
		elseif (isa(T.R,Missing)) && !(isa(T.L,Missing))
			push!(Q,T.L)
			push!(Q,Missing())
			print(T.O*Z)
			M=M*"|"*PP*"\\"*PPP
		elseif !(isa(T.R,Missing)) && (isa(T.L,Missing))	
			push!(Q,Missing())
			push!(Q,T.R)
			print(T.O*Z)
			M=M*"|"*PP*"\\"*PPP
		else
			print(T.O*Z)
			push!(Q,Missing());
			push!(Q,Missing());
			M=M*"|"*PP*"\\"*PPP
		end
	end
	println("");
	if i != 2
		println(M)
	end
	TreePrint(Q,i-1)
end
function Parse(S)
	println("Function: ",S);
	L = []
	for s in S
		append!(L,s)
	end
	println("Listed expression: ",L)
end
Parse(S)
TreeEx()
