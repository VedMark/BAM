#!/bin/bash

function associate_left() {
	echo; echo associating of the sample "$1"

	./BAM l "$1" "$2" out.txt && cat "$2"
}

function associate_right() {
	echo; echo associating of the sample "$1"

	./BAM r "$1" "$2" out.txt && cat "$2"
}

./BAM c 8 6 5 out.txt

./BAM a samples/A1.txt samples/B1.txt out.txt
./BAM a samples/A2.txt samples/B2.txt out.txt
./BAM a samples/A3.txt samples/B9.txt out.txt

associate_left samples/A1.txt A1.txt
associate_right samples/B1.txt B1.txt
associate_left samples/A2.txt A2.txt
associate_right samples/B2.txt B2.txt
associate_left samples/A9.txt A3.txt
associate_right samples/B9.txt B3.txt
