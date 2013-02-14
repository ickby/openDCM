
template< typename T >
void foo(T t) {};

extern template void foo<int> (int t) {};

int main()
{
    return 0;
}

