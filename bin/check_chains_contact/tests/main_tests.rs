use add;

mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }

    fn test_add() {
        let result = add(5, 6);
        assert_eq!(result, 11);
    }
}
