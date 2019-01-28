function search(input) {
  // Declare variables
  var elt, filter, uls, li, a, i;
  elt = input.parentElement;
  filter = input.value.toUpperCase();
  uls = elt.getElementsByClassName('thumbs');

  // Loop through all list items, and hide those who don't match the search query
  for (j = 0; j < uls.length; j++) {
    li = uls[j].getElementsByTagName('li');
    for (i = 0; i < li.length; i++) {
      a = li[i].getElementsByTagName("a")[0];
      if (a.innerHTML.toUpperCase().indexOf(filter) > -1) {
        li[i].style.display = "";
      } else {
        li[i].style.display = "none";
      }
    }
  }
}
