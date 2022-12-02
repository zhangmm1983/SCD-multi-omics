data = read.table()
col = colorRampPalette(c("navy", "white", "firebrick3"))(100)
tl.cex = 5
mode = "upper"
order = "hclust"
title = "\n Correlation"
tl.col = "black"
diag = F
sig.level = c("", "", "")
insig = "label_sig"
pch.cex = 6 
number.cex = 6 

data2 <- psych::corr.test(x = t(data), y = NULL, adjust = "none",ci = FALSE)

data2 <- corrplot::corrplot(
corr = data$r,
type = mode,
order = order,
title = title,
tl.col = tl.col,
tl.cex = tl.cex,
col = col,
diag = diag,
p.mat = data$p.adj,
insig = insig,
pch.cex = pch.cex,
number.cex = number.cex,
sig.level = sig.level)
