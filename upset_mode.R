
r = 1
h = 0.5 + 0.06
top_circle = function(from = 0, to = 2*pi) {
	theta = seq(from, to, length = 100)
	x = cos(theta)*r
	y = sin(theta)*r
	y = y + h
	data.frame(x, y)
}

left_circle = function(from = 0, to = 2*pi) {
	theta = seq(from, to, length = 100)
	x = cos(theta)*r
	y = sin(theta)*r
	x = x - h*sqrt(3)/2
	y = y - h/2
	data.frame(x, y)
}

right_circle = function(from = 0, to = 2*pi) {
	theta = seq(from, to, length = 100)
	x = cos(theta)*r
	y = sin(theta)*r
	x = x + h*sqrt(3)/2
	y = y - h/2
	data.frame(x, y)
}
library(grid)

grid.lines = function(df, ...) {
	grid::grid.lines(df[, 1], df[, 2], ..., default.units = "native")
}
grid.polygon = function(df, ...) {
	grid::grid.polygon(df[, 1], df[, 2], ..., default.units = "native")
}
make_venn = function(top, left, right, mode, fill = "red", col = "black") {
	if(missing(left) && missing(right)) {
		left = top[2]
		right = top[3]
		top = top[1]
	}
	x = c(-1.5, 0, 1.5)
	if(mode == "upset") {
		pushViewport(viewport(xscale = c(-2, 2), yscale = c(-2, 2)))
	} else {
		pushViewport(viewport(xscale = c(-2, 2), yscale = c(-2, 2),
			width = unit(1, "snpc"), height = unit(1, "snpc")))
	}
	if(top && left && right) { # 111
		if(mode == "distinct" || mode == "intersect") {
			grid.polygon(rbind(top_circle(-1/3*pi, -2/3*pi),
				               right_circle(pi, 2/3*pi),
				               left_circle(1/3*pi, 0)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "union") {
			grid.polygon(rbind(top_circle(pi, 0),
				               right_circle(1/3*pi, -2/3*pi),
				               left_circle(-1/3*pi, -4/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else {
			grid.points(x, rep(0, 3), pch = 16)
			grid::grid.lines(c(x[1], x[3]), c(0, 0), gp = gpar(lwd = 3), default.units = "native")
		}
	} else if(top && left && !right) { #110
		if(mode == "distinct") {
			grid.polygon(rbind(top_circle(pi, 4/3*pi),
				               right_circle(pi, 2/3*pi),
				               left_circle(1/3*pi, 2/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "intersect") {
			grid.polygon(rbind(top_circle(pi, 5/3*pi),
				               left_circle(0, 2/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "union") {
			grid.polygon(rbind(top_circle(-1/3*pi, pi),
				               left_circle(2/3*pi, 2*pi)),
			             gp = gpar(fill = fill, col = col))
		} else {
			grid.points(x, rep(0, 3), pch = 16, gp = gpar(col = c("black", "black", "#CCCCCC")))
			grid::grid.lines(c(x[1], x[2]), c(0, 0), gp = gpar(lwd = 3), default.units = "native")
		}
	} else if(top && !left && right) { #101
		if(mode == "distinct") {
			grid.polygon(rbind(top_circle(0, -1/3*pi),
				               left_circle(0, 1/3*pi),
				               right_circle(2/3*pi, 1/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "intersect") {
			grid.polygon(rbind(top_circle(0, -2/3*pi),
				               right_circle(pi, 1/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "union") {
			grid.polygon(rbind(top_circle(4/3*pi, 0),
				               right_circle(1/3*pi, -pi)),
			             gp = gpar(fill = fill, col = col))
		} else {
			grid.points(x, rep(0, 3), pch = 16, gp = gpar(col = c("black", "#CCCCCC", "black")))
			grid::grid.lines(c(x[1], x[3]), c(0, 0), gp = gpar(lwd = 3), default.units = "native")
		}
	} else if(!top && left && right) { #011
		if(mode == "distinct") {
			grid.polygon(rbind(top_circle(-1/3*pi, -2/3*pi),
				               right_circle(pi, 4/3*pi),
				               left_circle(-1/3*pi, 0)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "intersect") {
			grid.polygon(rbind(left_circle(1/3*pi, -1/3*pi),
				               right_circle(-2/3*pi, -4/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "union") {
			grid.polygon(rbind(left_circle(1/3*pi, 5/3*pi),
				               right_circle(-2/3*pi, 2/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else {
			grid.points(x, rep(0, 3), pch = 16, gp = gpar(col = c("#CCCCCC", "black", "black")))
			grid::grid.lines(c(x[2], x[3]), c(0, 0), gp = gpar(lwd = 3), default.units = "native")
		}
	} else if(top && !left && !right) { #100
		if(mode == "distinct") {
			grid.polygon(rbind(top_circle(0, pi),
				               left_circle(2/3*pi, 1/3*pi),
				               right_circle(2/3*pi, 1/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "intersect" || mode == "union") {
			grid.polygon(rbind(top_circle()),
			             gp = gpar(fill = fill, col = col))
		} else {
			grid.points(x, rep(0, 3), pch = 16, gp = gpar(col = c("black", "#CCCCCC", "#CCCCCC")))
		}
	} else if(!top && left && !right) { #010
		if(mode == "distinct") {
			grid.polygon(rbind(top_circle(pi, 4/3*pi),
				               right_circle(pi, 4/3*pi),
				               left_circle(-1/3*pi, -4/3*pi)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "intersect" || mode == "union") {
			grid.polygon(rbind(left_circle()),
			             gp = gpar(fill = fill, col = col))
		} else {
			grid.points(x, rep(0, 3), pch = 16, gp = gpar(col = c("#CCCCCC", "black", "#CCCCCC")))
		}
	} else if(!top && !left && right) { #001
		if(mode == "distinct") {
			grid.polygon(rbind(top_circle(-1/3*pi, 0),
				               right_circle(1/3*pi, -2/3*pi),
				               left_circle(-1/3*pi, 0)),
			             gp = gpar(fill = fill, col = col))
		} else if(mode == "intersect" || mode == "union") {
			grid.polygon(rbind(right_circle()),
			             gp = gpar(fill = fill, col = col))
		} else {
			grid.points(x, rep(0, 3), pch = 16, gp = gpar(col = c("#CCCCCC", "#CCCCCC", "black")))
		}
	}
	if(mode != "upset") {
		grid.lines(top_circle(), gp = gpar(col = col))
		grid.lines(left_circle(), gp = gpar(col = col))
		grid.lines(right_circle(), gp = gpar(col = col))
	}
	popViewport()
}

grid.newpage()
pushViewport(viewport(layout = grid.layout(nc = 4, nr = 8,
	height = unit.c(grobHeight(textGrob("A"))*2, unit(rep(1, 7), "null")))))
comb_mat = cbind(c(1, 1, 1),
	             c(1, 1, 0),
	             c(1, 0, 1),
	             c(0, 1, 1),
	             c(1, 0, 0),
	             c(0, 1, 0),
	             c(0, 0, 1))
all_modes = c("upset", "distinct", "intersect", "union")
for(i in 2:4) {
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
	grid.text(all_modes[i])
	popViewport()
}
for(i in 1:4) {
	for(j in 1:7) {
		pushViewport(viewport(layout.pos.row = j+1, layout.pos.col = i))
		make_venn(comb_mat[, j], mode = all_modes[i], fill = "black", col = "#808080")
		popViewport()
	}
}
popViewport()





