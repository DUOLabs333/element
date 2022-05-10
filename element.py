#!/usr/bin/env python3

# < include pyparsing/pyparsing.py >

from pyparsing import (Suppress, Word, nums, alphas, Regex, Forward, Group, 
						Optional, OneOrMore, ParseResults)


from collections import defaultdict

# < include "elements.json" pt >

# < include "table.txt" table >

table=table.decode()


import json
import readline
from pprint import pprint
import copy
import sys
import math
from decimal import Decimal
pt=json.loads(pt)

elements=list(map(str.lower, list(pt.keys())))

symbols=[pt[_]["Symbol"].lower() for _ in list(pt.keys())]

keepnames=["Symbol","Protons","Electrons","Neutrons", "Period", "Group", "Electronic configuration","Atomic mass", "Atomic number", "Atomic radius", "Density (g/cm^3)", "Specific heat (J/K)", "Electronegativity", "Melting point (K)", "Boiling point (K)", "First ionization (eV)"]

keepkeys=["Symbol","NumberofProtons","NumberofElectrons","NumberofNeutrons","Period","Group","ElectronicConfiguration", "AtomicMass","AtomicNumber","AtomicRadius","Density","SpecificHeat","Electronegativity","MeltingPoint","BoilingPoint","FirstIonization"]

def getResult(Element):
	
	if isinstance(Element,str):
		temp_dict={}
		for key in keepkeys:
			temp_dict[keepnames[keepkeys.index(key)]]=pt[Element][key]
			
		result={}
		result['Name']=Element
		result=result | temp_dict
	else:
		result={}
		result['Compound']=''.join([str(_[0])+str(_[1]) for _ in Element])
		result['Molar Mass']=float(sum([Decimal(_[1]*float(getResult(parseElement(_[0]))['Atomic mass'])) for _ in Element]))
	return result
	

def parseElement(Input):
	if Input.lower() in elements:
		return Input.title()
	elif Input.lower() in symbols:
		return elements[symbols.index(Input.lower())].title()
	else:
		raise Exception

def parseCompound(Input):
	LPAR,RPAR = map(Suppress,"()")
	integer = Word(nums)
	
	# add parse action to convert integers to ints, to support doing addition 
	# and multiplication at parse time
	integer.setParseAction(lambda dummy1,dummy2,t:int(t[0]))
	
	element = Word(alphas.upper(), alphas.lower())
	# or if you want to be more specific, use this Regex
	# element = Regex(r"A[cglmrstu]|B[aehikr]?|C[adeflmorsu]?|D[bsy]|E[rsu]|F[emr]?|"
	#                 "G[ade]|H[efgos]?|I[nr]?|Kr?|L[airu]|M[dgnot]|N[abdeiop]?|"
	#                 "Os?|P[abdmortu]?|R[abefghnu]|S[bcegimnr]?|T[abcehilm]|"
	#                 "Uu[bhopqst]|U|V|W|Xe|Yb?|Z[nr]")
	
	# forward declare 'formula' so it can be used in definition of 'term'
	formula = Forward()
	
	term = Group((element | Group(LPAR + formula + RPAR)("subgroup")) + 
					Optional(integer, default=1)("mult"))
	
	# define contents of a formula as one or more terms
	formula << OneOrMore(term)
	
	
	# add parse actions for parse-time processing
	
	# parse action to multiply out subgroups
	def multiplyContents(dummy1,dummy2,tokens):
		t = tokens[0]
		# if these tokens contain a subgroup, then use multiplier to
		# extend counts of all elements in the subgroup
		if t.subgroup:
			mult = t.mult
			for term in t.subgroup:
				term[1] *= mult
			return t.subgroup
	term.setParseAction(multiplyContents)
	
	# add parse action to sum up multiple references to the same element
	def sumByElement(w,x,tokens):
		elementsList = [t[0] for t in tokens]
	
		# construct set to see if there are duplicates
		duplicates = len(elementsList) > len(set(elementsList))
	
		# if there are duplicate element names, sum up by element and
		# return a new nested ParseResults
		if duplicates:
			ctr = defaultdict(int)
			for t in tokens:
				ctr[t[0]] += t[1]
			return ParseResults([ParseResults([k,v]) for k,v in ctr.items()])
	formula.setParseAction(sumByElement)
	return formula.parseString(Input)

def parseInput(Input):
	parseCompound(Input)
	try:
		if Input.islower():
			result=parseElement(Input)
		else:
			result=parseCompound(Input)
	except:
		print("Incorrect Element/Molecule")
		return
	pprint(getResult(result),sort_dicts=False,width=1)
	
def repl():
	while True:
		Input=input('>> ')
		parseInput(Input)

if len(sys.argv)==1:
	try:
	   repl()
	except (EOFError,KeyboardInterrupt):
		print()
		exit()
else:
	if sys.argv[-1]=="--table":
		print(table)
		exit()
	Input=sys.argv[1]
	parseInput(Input)

	