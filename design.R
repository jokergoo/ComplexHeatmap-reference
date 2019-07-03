

single_heatmap_layout = function() {
grid.newpage()

th = convertHeight(grobHeight(textGrob("a")) * 3, "mm")
pushViewport(viewport(layout = grid.layout(nr = 3, nc = 3,
	width = unit.c(th*4, unit(1, "null"), th*4),
	height = unit.c(th*4, unit(1, "null"), th*4)),
	gp = gpar(fontsize = 12), 
	width = unit(1, "npc") - unit(4, "mm"),
	height = unit(1, "npc") - unit(4, "mm")))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
name = c("column annotations", "column names", "dendrogram", "title")
for(i in 1:4) {
	pushViewport(viewport(y = i/4, height = 1/4, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
name = rev(c("column annotations", "column names", "dendrogram", "title"))
for(i in 1:4) {
	pushViewport(viewport(y = i/4, height = 1/4, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
name = rev(c("row annotations", "row names", "dendrogram", "title"))
for(i in 1:4) {
	pushViewport(viewport(x = i/4, width = 1/4, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
name = c("row annotations", "row names", "dendrogram", "title")
for(i in 1:4) {
	pushViewport(viewport(x = i/4, width = 1/4, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
name = rbind(c("matrix, row slice 1\ncolumn slice1", "matrix, row slice 1\ncolumn slice2"), 
	         c("matrix, row slice 2\ncolumn slice1", "matrix, row slice 2\ncolumn slice2"))
for(i in 1:2) {
	for(j in 1:2) {
		pushViewport(viewport(x = i/2, y = j/2, width = 1/2, height = 1/2, just = c("right", "top")))
		pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"),
			height = unit(1, "npc") - unit(4, "mm")))
		grid.rect(gp = gpar(col = "red", fill = "#FF000020"))
		grid.text(name[i, j])
		popViewport()
		popViewport()
	}
}
popViewport()

popViewport()

}




heatmap_list_layout = function(direction = "horizontal") {
grid.newpage()

th = convertHeight(grobHeight(textGrob("a")) * 3, "mm")
pushViewport(viewport(layout = grid.layout(nr = 3, nc = 3,
	width = unit.c(th*2, unit(1, "null"), th*2),
	height = unit.c(th*2, unit(1, "null"), th*2)),
	gp = gpar(fontsize = 12),
	width = unit(1, "npc") - unit(4, "mm"),
	height = unit(1, "npc") - unit(4, "mm")))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
name = c("title", "legends")
for(i in 1:2) {
	pushViewport(viewport(y = i/2, height = 1/2, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
name = rev(c("title", "legends"))
for(i in 1:2) {
	pushViewport(viewport(y = i/2, height = 1/2, just = "top"))
	grid.rect()
	grid.text(name[i])
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
name = rev(c("title", "legends"))
for(i in 1:2) {
	pushViewport(viewport(x = i/2, width = 1/2, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
name = c("title", "legends")
for(i in 1:2) {
	pushViewport(viewport(x = i/2, width = 1/2, just = "right"))
	grid.rect()
	grid.text(name[i], rot = 90)
	popViewport()
}
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
name = c("heatmap1", "heatmap2", "row annotation")
for(i in 1:3) {
	if(direction == "horizontal") {
		name = c("heatmap1", "heatmap2", "row annotation")

		pushViewport(viewport(x = i/3, width = 1/3, just = "right"))
		pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"),
			height = unit(1, "npc") - unit(4, "mm")))
		grid.rect(gp = gpar(col = "green", fill = "#00FF0020"))
		grid.text(name[i], rot = 90)
		popViewport()
		popViewport()
	} else {
		name = rev(c("heatmap1", "heatmap2", "column annotation"))

		pushViewport(viewport(y = i/3, height = 1/3, just = "top"))
		pushViewport(viewport(width = unit(1, "npc") - unit(4, "mm"),
			height = unit(1, "npc") - unit(4, "mm")))
		grid.rect(gp = gpar(col = "red", fill = "#FF000020"))
		grid.text(name[i])
		popViewport()
		popViewport()
	}
}
popViewport()

popViewport()

}







