<?php

define('STORE', 'ratings.dat');

function put_contents($file,$content) {
	$f = fopen($file,"w");
	fwrite($f,$content);
	fclose($f);
}

if(isset($_REQUEST["rating"])) {
	$rating = $_REQUEST["rating"];
	$storedRatings = unserialize(file_get_contents(STORE));
	$storedRatings[] = $rating;
	put_contents(STORE, serialize($storedRatings));
	$average = round(array_sum($storedRatings) / count($storedRatings), 2);
	$count = count($storedRatings);
	$xml = "<ratings><average>$average</average><count>$count</count></ratings>";
	header('Content-type: text/xml');
	echo $xml;
}

?>