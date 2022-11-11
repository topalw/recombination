# A few thoughts on smc++ in CH  

Sequential markovian coalescent methods tend to show a similar pattern of Ne oscillation in the past for many species and are widely 
considered untrustworthy at least unofficially!
SMC++ has been shown to perform well in simulations (ex. stdpopsim paper) so I will assume that the results it outputs are (at least partially) reliable. 
We will see later that they are actually not so much - at least by a factor of 2-3 - which in my humble opinion is not that bad. 

So assuming that the curve we see in CH.pdf is the actual history of the population we can assume that the model would identify the out of refugium 
bottleneck (oor-bneck) of the population (or the glacial age contraction). Care that we do not necessarily expect a population contraction at the onset of the 
glaciation since the new colonizers of the population might be unrelated to the owls that experienced the onset of the event (which might have been extinct). 
However we do expect a population contraction at the colonizing stage unless the migration from the source is constant and big (which is not impossible) or 
total population replacement happened in the time since (which is a bit more sketchy). 
But assume that this oor-bneck is what we see in the curve. The temporal location of it is at 2e3 to 4e3 generations in the past. This translates to 
5-10'000 years ago which actually coincides well with the oor-bneck or even the glacial expansion at the onset of the glaciation. 
Now there are a few caveats associated withw SMC or any such method. 
One could be the input unrelatedness and how well you sample to pop coalescence tree (with 93 inds unrelated at 0.03 assume we do).

Another is the neccessity of a mutation rate which in our case comes from the flycatcher family paper which is a bit sketchy. 
The mutation rate essentialy scales the diversity present in the haplotypes to an estimation of Ne through the 4Nem parameter. 
If that value is correct we expect the Ne used in pyrho to be accurate (or at least in a ballpark of reality). 
What happens in pyrho is that almost consistently, both for GR and CH we see an underestimation of recombination by a factor of 2. 
This number comes from comparing the cM length from LepMap to the cM length of pyrho. 
Pyrho uses Ne from smcpp to scale the population parameter rho = 4Ne*theta* to *theta*. 
Essentially this means that if the LM3 estimates of cM are accurate (they better be!) then 4Ne has been overestimated by 2. 
If we look back at the curve we see a contemporary (or actually a bit old (1k gens!)) Ne in CH of 100'000 which is a bit absurd. 
Even the half of it - post pyrho Ne correction is a bit crazy but at least half as bad. 
The thought at the start of this rant is that probably mutation rate that scales the Ne y axis of smc++ is off and that leads to erroneous 
estimates of Ne for all our pops. In fact with the CH example the mutation rate is underestimated (lower mu means less need of Ne to increase pi). 
However the Ne estimation of GR is much more reasonable ~3e4 (still large) and the question is what prompts such an increased Ne in CH 
(could the individual relatedness do sth like this?).
