./bpls.sh | grep contact_map | tr -s [:space:] | cut -d' ' -f4 | cut -d'*' -f1 | perl -e '@a=<>;use List::Util qw(sum);print sum(@a)'
echo ""
