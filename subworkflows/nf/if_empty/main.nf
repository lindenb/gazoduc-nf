workflow IF_EMPTY {
take:
	channel_A
	channel_B
main:
	fa = channel_A.map{[1,it]}
	fb = channel_B.map{[2,it]}

	// min value
	m = fa.first()
		.mix(fb.first())
		.map{it[0]}
		.min()

	fc = fa.mix(fb)
		.combine(m)
		.filter{idx,row,min_v->idx==min_v}
		.map{idx,row,min_v->row}
emit:
	output = fc
}
