library(motifRG)
data(ctcf.seq)
data(control.seq)
### concatenate the foreground, background sequences
all.seq <- append(ctcf.seq, control.seq)
### specify which sequences are foreground, background. 
category <- c(rep(1, length(ctcf.seq)), rep(0, length(control.seq))) 

### find motifs
ctcf.motifs <- findMotif(all.seq=all.seq, category=category, max.motif=3)

###print the summary of motifs
summaryMotif(ctcf.motifs$motifs, ctcf.motifs$category)

###print sequence logo of the first motif
plotMotif(ctcf.motifs$motifs[[1]]@match$pattern)

###Create a table for motifs for a latex document
motifLatexTable(main="CTCF motifs", ctcf.motifs)

###Find a refined PWM model given the motif matches as seed
pwm.match <- refinePWMMotif(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
library(seqLogo)
seqLogo(pwm.match$model$prob)


### Motifs found by findMotif tend to be relatively short, as longer and more specific
### motif models do not necessarily provide better discrimination of foreground background 
### vs background if they are already well separated.
### However, one can refine and extend a PWM model given the motif matches by findMotif as seed
### for more specific model. 
pwm.match.extend <- refinePWMMotifExtend(ctcf.motifs$motifs[[1]]@match$pattern, ctcf.seq)
seqLogo(pwm.match.extend$model$prob)
###plot the dinucleotide logo fo PWM match
plotMotif(pwm.match.extend$match$pattern)

###Show dependency of adjacent positions, green for enriched pair, red for depleted pair
plotMotif(pwm.match.extend$match$pattern, logodds=T)







