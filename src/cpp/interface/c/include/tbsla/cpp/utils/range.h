#ifndef TBSLA_CINTERFACE_CPP_range
#define TBSLA_CINTERFACE_CPP_range

#ifdef __cplusplus
extern "C" {
#endif

static long long int lnv(long long int s, long long int l, long long int g)
{
  long long int n = s / g;
  long long int mod = s % g;
  if (l < mod)
    n++;
  return n;
}

static long long int pflv(long long int s, long long int l, long long int g)
{
  long long int mod = s % g;
  long long int n = lnv(s, l, g) * l;
  if (l >= mod)
    n += mod;
  return n;
}

#ifdef __cplusplus
}
#endif

#endif /* TBSLA_CINTERFACE_CPP_Vector */
